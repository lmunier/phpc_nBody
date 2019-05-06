//
// Created by lmunier on 06/05/19.
//

#ifndef PROJECT_PARTICLE_HPP
#define PROJECT_PARTICLE_HPP

#include "constants.hpp"
#include "Vector.hpp"

template<typename Type>
class Particle {
public:
    Particle (float mass, Type pos) : pos(pos), m(mass) {}

    // Return value of updated to know if this particle is already updated to new state
    bool is_updated() {return this->updated;}

    // Return private variable pos
    Type get_pos() {return this->pos;}

    // Return private variable vel
    Type get_vel() {return this->vel;}

    // Return private variable vel
    Type get_load() {return this->load;}


private:
    // Set boolean to keep a state if updated already or not
    bool updated = false;

    // Store position, velocity and load of the particle in space domain
    Type pos;
    Type vel;
    Type load;

    // Store mass of the particle
    float m = 0.0f;
};


#endif //PROJECT_PARTICLE_HPP
