/**
 * @file Particle.hpp
 * @author Munier Louis
 * @date 03.06.19
 * @version 3.0
 *
 * Particle file to implement a class named particle,
 */

#ifndef PHPC_NBODY_PARTICLE_HPP
#define PHPC_NBODY_PARTICLE_HPP

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 * @include Vector.hpp custom library to have minimal vector implementation
 * @include Tree.hpp custom class to have mother class link
 */
#include "constants.hpp"

/**
 * Class Particle, to store all the information of a given particle.
 *
 * @class Particle
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 */
struct Particle {
    float _pos[NB_DIM] = {0.0f};                    /**< position of the given particle */
    float _vel[NB_DIM] = {0.0f};                    /**< velocity of the given particle */
    float _load[NB_DIM] = {0.0f};                   /**< load on the given particle */
    float _m = 0.0f;                                /**< mass of the given particle */
};

/**
 * Update velocity and position of a given particle for a given load on it. Reset load after update.
 */
// void update_vel_pos() {
//     float new_velocity[NB_DIM] = this->get(LOAD) * (DELTA_T / this->get_mass()) + this->get(VEL);
//     float new_position[NB_DIM] = new_velocity * DELTA_T + this->get(POS);

//     this->set(VEL, new_velocity);
//     this->set(POS, new_position);
//     this->set(LOAD, Type());
// }

/**
 * Test if a given particle is out of the overall area and return true if it is the case, false otherwise.
 *
 * @return true if particle out of its boundaries, false otherwise
 */
// bool is_out_boundaries(Type dim) {
//     Type pos = this->get(POS);

//     if (abs(pos.x) > 0.5 * dim.x)
//         return true;

//     if (abs(pos.y) > 0.5 * dim.y)
//         return true;

// #if NB_DIM == DIM_3
//     if (abs(pos.z) > 0.5 * dim.z)
//         return true;
// #endif
//     return false;
// }

#endif //PHPC_NBODY_PARTICLE_HPP
