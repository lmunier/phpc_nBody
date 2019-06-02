/**
 * @file constants.hpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
 *
 * Constants file to have all the desired constants of the overall project.
 */

#ifndef PROJECT_CONSTANTS_HPP
#define PROJECT_CONSTANTS_HPP

/**
 * Different includes to have all the c++ functions/library needed.
 */
#include <iostream>            /**< input/output stream objects */
#include <fstream>             /**< input/output stream class to operate on files */
#include <string>              /**< string library manipulation */
#include <cstdlib>             /**< C standard general utilities library */
#include <cstdint>             /**< integer types */
#include <random>              /**< introduces random number generation facilities */
#include <vector>              /**< dynamic arrays */
#include <math.h>              /**< different maths functions */

#include <chrono>              /**< deal with time */

/**
 * @namespace std::chrono to simplify notation on time duration computing
 */
using namespace std::chrono;

/**
 * Different defines to have constants for the overall project.
 */
#define PRINT                          /**< if it is defined, print value in terminal and csv file */

#define G 6.67408e-11f                  /**< gravitational general constant */
#define EPSILON 0.000005f               /**< to avoid collision between particles */
#define BH_THETA 0.5f                   /**< Barnes-Hut threshold to consider particles as a unique one */

#define DIM_2 2                        /**< number of dimension in 2D */
#define DIM_3 3                        /**< number of dimension in 3D */
#define NB_DIM DIM_3                   /**< number of dimension chosen for the project */

#define NB_PARTICLES 1000                /**< number of particles for the project */
#define SIDE 1000                      /**< side of the area considered for the project */
#define SHIFT SIDE/3.0f                 /**< shift each particles to unbalanced probability of particles position */
#define OCCUPATION_PERC 0.5f            /**< percentage of occupation to avoid particle go easily outside boundaries */
#define MASS_MAX 10e10f                 /**< maximum of mass value for a particle */

#define DELTA_T 0.01f                   /**< duration of each update timestep */
#define ITERATIONS 1000                /**< number of iterations to solve nBody problem */

/**
 * @enum property of the particle
 */
enum property {
    POS,          /**< enum value for position property */
    VEL,          /**< enum value for velocity property */
    LOAD,         /**< enum value for load property */
};

#endif //PROJECT_CONSTANTS_HPP
