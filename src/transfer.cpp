#include "transfer.h"

Transfer::Transfer(){

}

Transfer::~Transfer(){

}

void Transfer::add_planets(Planet* p1, Planet* p2){
    /**
     * @brief Adds the two reference plantes
     */
    //std::cout << "Added  1:" << p1 << " Added  2: "<< p2 << std::endl;
    this->p1 = p1;
    this->p2 = p2;
}

void Transfer::compute_transfer(float T){
    /**
     * @brief Computes the transfer between both planets in the class.
     */
    orbit::lambert(this->p1->r_eph, this->p2->r_eph, T, MU_SUN, this->v_dep, this->v_arr);
}

void Transfer::print() const{

}
