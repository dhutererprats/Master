#include "genetic.h"

ProblemDefinition::ProblemDefinition(const float _dep) : departure(_dep){
}

ProblemDefinition::~ProblemDefinition(){
}

void ProblemDefinition::add_planet(int _p, float min, float max){
    
    this->planets.push_back(Planet(_p));
    if(min >= max){
        throw timeRange();
    }

    this->timeWindows.push_back(std::make_pair(min, max));
}

void ProblemDefinition::set_max_time(float t){
    this->time_max = t;
}


/*
====================================================================================================================
                                                I N D I V I D U A L                                               ==
====================================================================================================================
*/

Individual::Individual(ProblemDefinition* prob){
    this->problem = prob;
    this->cost = 0.f;
    this->fitness = 0.f;
    this->totalDV = 0.f;
}

Individual::~Individual(){

}

std::string Individual::getChromosome() const{
    /**
     * @brief Get full chromosome in a bitstring
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
    if(at >= this->flyTimes.size() || at < 0){
        throw "Gene position out of range";
    }
    //std::cout << this->flyTimes.at(at) << std::endl;
    return time2bitStr(this->flyTimes.at(at));
}

void Individual::setChromosome(std::string chromo){
    /**
     * @brief Given a full chromosome, convert to flight times and update it.
     */
    for(unsigned int i = 0; i < this->flyTimes.size(); i++){
        this->setGene(chromo.substr(i * (MAX_BIT_SIZE + MAX_FRACTIONAL_BIT_SIZE),(MAX_BIT_SIZE + MAX_FRACTIONAL_BIT_SIZE)), i);
    }   
}

void Individual::setGene(std::string gene, int at){
    /**
     * @brief Given a gene and position convert to flytime and update.
     */
    if(at >= this->flyTimes.size() || at < 0){
        throw "Gene position out of range";
    }
    float new_time = bitStr2Time(gene);
    new_time = (new_time >= this->problem->timeWindows.at(at).first) ? new_time: this->problem->timeWindows.at(at).first;
    new_time = (new_time <= this->problem->timeWindows.at(at).second) ? new_time : this->problem->timeWindows.at(at).second;
    this->flyTimes.at(at) = new_time;
}


void Individual::mate(const Individual& partner, int crossType){
    /**
     * @brief Matting function that creates a new individual. 
     * The new individual is actually the current one, value on him are update. Parent -> Transforms to "child".
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
        const int r_gene = rand() % this->flyTimes.size();
        std::string p1_g = time2bitStr(this->flyTimes.at(r_gene));
        std::string p2_g = time2bitStr(partner.flyTimes.at(r_gene));
        std::string child_g = uniformBitstrCross(p1_g, p2_g);

        this->setGene(child_g, r_gene);
    }
    else if(crossType == CROSS_SINGLE_POINT){
        // From the full chromosome. At a given random point, interchange the parent 1 bits by the partner's bits.
        std::string p1_chromo = this->getChromosome();
        std::string p2_chromo = partner.getChromosome();
        
        const int cut_at = std::rand() % p1_chromo.size();
        std::string child;
        child.append(p1_chromo.substr(0, cut_at)); // Start to cut_at.
        child.append(p2_chromo.substr(cut_at)); // from cut_at pos to end.

        this->setChromosome(child);
    }
    else if(crossType == CROSS_DOUBLE_POINT){
        // From the full chromosome. At a given random point 1, interchange the parent 1 bits by the partner's bits. Do the same a second time at point 2.
        std::string p1_chromo = this->getChromosome();
        std::string p2_chromo = partner.getChromosome();
        
        const int cut_at_1 = rand_rng(1, p1_chromo.size() - 2); // From above 1, to at least leave 1 for the second cut.
        const int cut_at_2 = rand_rng(cut_at_1 + 1,  p1_chromo.size() - 1); // From at least one bit after cut one to end leave at least one for this cut.
        std::string child;
        child.append(p1_chromo.substr(0, cut_at_1)); 
        child.append(p2_chromo.substr(cut_at_1, cut_at_2 - cut_at_1));
        child.append(p1_chromo.substr(cut_at_2));

        this->setChromosome(child);
    }
    else if (crossType == CROSS_PERSONALIZED){
        // Times are crossed staticly. As Half points (centroids) between parents depending on their fitness.
        // Done per gene to have different solution, otherwise per 2 parents there exist only one solution = 1 child.
        const int r_gene = rand() % this->flyTimes.size();
        float fit_w = (this->fitness / (this->fitness + partner.fitness)); 
        int sign = (this->flyTimes.at(r_gene) > partner.flyTimes.at(r_gene))? -1 : 1;
        float centroid_t = this->flyTimes.at(r_gene) + sign * std::fabs(this->flyTimes.at(r_gene) - partner.flyTimes.at(r_gene))*fit_w; 
        std::string chlild_t = time2bitStr(centroid_t);
        this->setGene(chlild_t, r_gene);

    }
    else{
        throw "Unkown crossover type!";
    }
    
}

void Individual::createMutation(){
    /**
     * @brief Generates mutation in the individual.
     */
    std::string chromo = this->getChromosome();
    int mutate_at = rand_rng(0, chromo.size() - 1);
    chromo.at(mutate_at) = (chromo.at(mutate_at) == '0')? '1' : '0';

    this->setChromosome(chromo);
}


void Individual::init(){
    /**
     * @brief Used to initizilize the individual to a random vector of times (within the windows).
     */
    for(const auto& window: this->problem->timeWindows){
        float t = window.first + rand() % (int)(window.second - window.first - 1);
        t += rand_d(); // Add decimals to inrease resolution.
        this->flyTimes.push_back(t);
    }
}


int Individual::getFlyTime() const{

}

void Individual::evaluate(){
    /**
     * @brief Solves the mga problem. Updates the cost.
     * 
     */

    // Per planets
    for(unsigned int i = 0; i < this->problem->planets.size() - 2; i++){
        double r1[3], v1[3], r2[3], v2[3], r3[3], v3[3];
        double T1 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + i + 1, 0.0); // Times are in format [Departure T, +ΔT1, +ΔT2,..]
        double T2 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + (i + 1) + 1, 0.0);
        double T3 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + (i + 2) + 1, 0.0);
        orbit::ephemeris(this->problem->planets.at(i).prm, T1, r1, v1);
        orbit::ephemeris(this->problem->planets.at(i + 1).prm, T2, r2, v2);
        orbit::ephemeris(this->problem->planets.at(i + 2).prm, T3, r3, v3);

        double v_dep1[3], v_arr2[3], v_dep2[3], v_arr3[3];

        orbit::lambert(r1, r2, (T2 - T1) * DAY2SEC, MU_SUN, v_dep1, v_arr2);
        orbit::lambert(r2, r3, (T3 - T2) * DAY2SEC, MU_SUN, v_dep2, v_arr3);

        double dV, delta, peri;
        orbit::patched_conic(v_arr2, v_dep2, v2, this->problem->planets.at(i + 1).mu, dV, delta, peri);

        // Cost update if departure too.
        // std::cout << "------------------------------------------" << std::endl;
        if(i == 0){
            // Departure is speed must be considered at infinity, accouting for the speed of earth, automaticlly given to the prove.
            double v_dep_at_infi[3];
            minus2(v_dep1, v1, v_dep_at_infi);
            this->updateDepartureCost(norm(v_dep_at_infi));
            this->totalDV += norm(v_dep_at_infi);
        }

        // Cost update
        this->updateCost(this->problem->planets.at(i + 1), dV, delta, peri, norm(v_arr2));
        this->totalDV += dV;
    }


    // Time penalty cost if the time is to big.
    float total_time = std::accumulate(this->flyTimes.begin(), this->flyTimes.end(), 0.0);
    if( total_time > this->problem->time_max){
        this->cost += K_TIME_PENALTY * (total_time - this->problem->time_max); 
    }
    this->fitness = this->cost;
}


void Individual::updateCost(const Planet& planet, double dV, double delta, double peri, double vin){
    /**
     * @brief Updates the current cost.
     * - Accumulates the delta V of the flyby.
     * - If the flyby is not feasable (too small peri for instance) it increases the cost a lot (to avoid solution evolving).
     */
    // TODO: change to also consider the perigee and maybe the angle too???
    this->cost += dV;

    // Penalty function.
    //this->cost += -2*std::log10(peri / (1.1 * planet.rad)); // Penalize low perigee radius. //FIXME: This is BS
    if(peri < 1.1 * planet.rad){
        this->cost += 1e10;
    }

    double rsoi = std::pow((planet.mass / MASS_SUN), 2/5) * planet.sun_dist;
    double E = ((vin * vin) / 2) - planet.mu/rsoi; // Penalize low velocities flybys (can lead to spacecraft planet capture).
    if(E < 0){
        this->cost += 1/vin;
    }
}

void Individual::updateDepartureCost(double dV){
    /**
     * @brief Adds departure cost.
     */
    this->cost += dV;
}




/*
====================================================================================================================
                                                P O P U L A T I O N                                               ==
====================================================================================================================
*/
 
Population::Population(const GenOperators params, ProblemDefinition* problem) : geParameters(params), 
                                                                                population(N_POPULATION, Individual(problem)), 
                                                                                newPopulation(N_POPULATION, Individual(problem)){

}

Population::~Population(){

}

void Population::inception(){
    /**
     * @brief Initializes each individual randomly.
     */
    for(auto& ind: this->population){
        ind.init();
        ind.evaluate();
    }
}

void Population::sortPopulation(){
    /**
     * @brief Sorts the population vector from fitest to less fit.
     */
    std::sort(this->population.begin(), this->population.end());
}

void Population::plotFitnessEvolution(){
    /**
     * @brief Plots the fitness evolution.
     */
    // TODO: Finsih me
    for(const auto& f: this->fitnessEvolution){
        std::cout << f << ",";
    }
    std::cout<<std::endl;
}

// ----- Genetic operators. -----
void Population::selection(){
    /**
     * @brief Selects the individuals that will undergo reporduction/crossover.
     */
    if(this->geParameters.selectionType == SELECTION_ROULETTE){
        float sum_adjusted_ft = 0.f;
        float roulette[N_POPULATION];
        
        for(unsigned int i = 0; i < N_POPULATION; i++){
            sum_adjusted_ft += 1/(1 + this->population.at(i).fitness);
        }
        // Roulette wheel
        for(unsigned int j = 0; j < N_POPULATION; j++){
            roulette[j] = (1/(1 + this->population.at(j).fitness)) / (sum_adjusted_ft);
        }
        // Do the selection
        float mean = 0;
        for(unsigned int k = 0; k < N_POPULATION - this->geParameters.elitism_n; k++){
            float r = rand_d();
            float curr = roulette[0];
            int idx = 0;
            while (r >= curr){
                idx++;
                curr += roulette[idx]; 
            }
            idx = (idx <= N_POPULATION - 1)? idx: N_POPULATION - 1; //Sanity-check (the >= for floating points is not excat enough! In the end they value the same.)
            this->newPopulation.push_back(this->population.at(idx));
        }        
    }
    else if(this->geParameters.selectionType == SELECTION_TOURNAMENT){
        // Tournament method
        float mean = 0.f;
        for(unsigned int i = 0; i < N_POPULATION - this->geParameters.elitism_n; i++){

            int rnd_idx[TOURNAMENT_N];
            for(int r = 0; r < TOURNAMENT_N; r++){
                rnd_idx[r] = rand_rng(0, N_POPULATION - 1);
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
    if(this->newPopulation.size() != N_POPULATION){
        throw "The population expected to be crossed is missing individuales. Size is different from N_POPULATION";
    }
    
    // From elitism size upwards, there are the selected individuals.
    for(unsigned int i = this->geParameters.elitism_n; i < N_POPULATION - 1; i = i+2){
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
     * @brief Mutate the generation.
     */
    for(unsigned int i = this->geParameters.elitism_n - 1; i < N_POPULATION -  this->geParameters.elitism_n; i++){
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
     * @brief Evaluates each individual
     */
    // for(auto& ind: this->population){
    //     ind.cost = 0;
    //     ind.fitness = 0;
    //     ind.totalDV = 0;
    //     ind.evaluate();
    //     ind.fitness = ind.cost; // delerte me, do somewhere else (TODO:);
    // }
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
    

}

void Population::evolveNewGenerationThreaded(int indx_start, int indx_end){
    /**
     * @brief Base evolution for multithreaded evolution
     */
    for(unsigned int i = indx_start; i < indx_end; i++){
        this->population.at(i).cost = 0;
        this->population.at(i).fitness = 0;
        this->population.at(i).totalDV = 0;
        this->population.at(i).evaluate();
        this->population.at(i).fitness = this->population.at(i).cost; // delerte me, do somewhere else (TODO:);
    }
}

void Population::runGeneration(){
    /**
     * @brief Runs the full genetic algorithm
     */

    int g = 0;
    this->newPopulation.clear();

    while (g < GEN_LIMIT){
        this->sortPopulation();
        // Record evolution of the fitness
        this->fitnessEvolution.push_back(this->population.at(0).fitness);

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


void Population::mateIndividuals(Individual parent1, Individual parent2){
    /**
     * @brief Mate two indivuals to produce to childs at the same time.
     * Child will be the oposite.
     */

}