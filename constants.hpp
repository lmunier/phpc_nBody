/**
 * @file constants.hpp
 * @author Munier Louis
 * @date 27.05.19
 * @version 1.0
 *
 * Constants file to have all the desired constants of the overall project.
 */

#ifndef PROJECT_CONSTANTS_HPP
#define PROJECT_CONSTANTS_HPP
/**
 * Different includes to have all the c++ functions/library needed.
 *
 * @namespace use std::chrono to simplify notation on time duration computing
 */
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <vector>
#include <math.h>

#include <chrono>
using namespace std::chrono;

/**
 * Different defines to have constants for the overall project.
 */
#define PRINT                  /**< if it is defined, print value in terminal and csv file to verify implementation */

#define G 6.67408              /**< gravitational general constant */
#define EPSILON 0.005          /**< to avoid collision between particles */

#define DIM_2 2                /**< number of dimension in 2D */
#define DIM_3 3                /**< number of dimension in 3D */
#define NB_DIM DIM_3           /**< number of dimension chosen for the project */

#define NB_PARTICLES 10        /**< number of particles for the project */
#define SIDE 100             /**< side of the area considered for the project */
#define SHIFT SIDE/3.0         /**< shift each particles position to unbalanced probability of particles position */
#define MASS_MAX 10            /**< maximum of mass value for a particle */

#define DELTA_T 0.01           /**< duration of each update timestep */
#define ITERATIONS 1000        /**< number of iterations to solve nBody problem */

/**
 * @enum property of the particle
 */
enum property {
    POS,          /**< enum value for position property */
    VEL,          /**< enum value for velocity property */
    LOAD,         /**< enum value for load property */
    MASS,         /**< enum value for mass property */
};

#endif //PROJECT_CONSTANTS_HPP
