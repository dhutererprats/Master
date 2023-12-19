#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>
#include <string>
#include <ctime>
#include <bitset>
#include <iostream>
#include "nlohmann/json.hpp"
#include "astro.h"


#define MAX_BIT_SIZE 11             // 2048 (max time per leg ~ 11 years).
#define MAX_FRACTIONAL_BIT_SIZE 3   //8 :: 1/8 = 0.125 of resolution for times
// Time value = [MAX_BIT_SIZE],[MAX_DECIMAL_BITS] 

void minus2(const double* v1, const double* v2, double* out);

double norm(const double* v1);

double distance(const double* v1, const double* v2);

double dot_prod(const double* v1, const double* v2);

double vec2(const double* v1);

void rotate_eph(double w, double W, double i, double P, double Q, double* vec);

void cross_prod(const double* v1, const double* v2, double* out);

double tofabn(const double &sigma,const double &alfa,const double &beta);

void vers(const double* v1, double* v2);

double x2tof(const double &x, const double &s, const double &c, bool prograde);

double deg2rad(const double x);

double rad2deg(const double x);

float rand_d();

int rand_rng(int min, int max);

std::string time2bitStr(const float x);

float bitStr2Time(const std::string x);

std::string uniformBitstrCross(const std::string s1, const std::string s2);


// Pagmo
void Conversion(const double*, double*, double*, const double &);

double Mean2Eccentric (const double &, const double &);

#endif //UTILITIES_H
