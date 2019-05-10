//
// Created by lmunier on 06/05/19.
//

#ifndef PROJECT_PARTICLE_HPP
#define PROJECT_PARTICLE_HPP

#include "constants.hpp"
#include "Vector.hpp"

#include <iostream>
#include <string>

template<typename Type>
class Particle {
public:
    explicit Particle (float mass = 0.0f, Type pos = Type()) : pos(pos), m(mass) {}

    // Return value of updated to know if this particle is already updated to new state
    bool is_updated() {return this->updated;}

    // Set update state fo the particle
    void update(bool state) {this->updated = state;}

    // Return private variable
    Type get(property p) {
        switch (p) {
            case POS :
                return this->pos;
            case VEL :
                return this->vel;
            case LOAD :
                return this->load;
            default:
                break;
        }
    }

    // Set private variable from particle
    Type set(property p, Type vec) {
        switch (p) {
            case POS :
                return this->pos = vec;
            case VEL :
                return this->vel = vec;
            case LOAD :
                return this->load = vec;
            default:
                break;
        }
    }


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
