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
#include <math.h>              /**< different maths functions */

#include <curand.h>

/**
 * Different defines to have constants for the overall project.
 */
//#define PRINT                          /**< if it is defined, print value in terminal and csv file */

#define G 6.67408e-11f                  /**< gravitational general constant */
#define EPSILON 0.005f               /**< to avoid collision between particles */

#define DIM_3 3                        /**< number of dimension in 3D */
#define NB_PARTICLES 500000                /**< number of particles for the project */
#define SIDE NB_PARTICLES                     /**< side of the area considered for the project */
#define SIDE_2f SIDE / 2.0f                  /**< 1/2 side of the area considered for the project */
#define SIDE_4f SIDE / 4.0f                   /**< 1/4 side of the area considered for the project */
#define SHIFT SIDE / 3.0f                     /**< Shift value to unbalanced uniformity of particles creation */
#define OCCUPATION_PERC 0.5f            /**< percentage of occupation to avoid particle go easily outside boundaries */
#define MASS_MAX 10e10f                 /**< maximum of mass value for a particle */

#define MAX_THREADS 1024

#define DELTA_T 0.01f                   /**< duration of each update timestep */
#define ITERATIONS 500                /**< number of iterations to solve nBody problem */

#endif //PROJECT_CONSTANTS_HPP
