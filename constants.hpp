//
// Created by lmunier on 06/05/19.
//

#ifndef PROJECT_CONSTANTS_HPP
#define PROJECT_CONSTANTS_HPP

#include <math.h>

#define PRINT

#define G 6.67408
#define epsilon 0.005

#define DIM_2 2
#define DIM_3 3
#define NB_DIM DIM_2

#define NB_PARTICLES 5000
#define SIDE 1000000
#define SHIFT SIDE/3.0
#define MASS_MAX 10

enum property {POS, VEL, LOAD, MASS};

#endif //PROJECT_CONSTANTS_HPP
