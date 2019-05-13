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
    explicit Particle (float mass = 0.0f, Type pos = Type()) : _pos(pos), _m(mass) {}

    // Return value of updated to know if this particle is already updated to new state
    bool is_updated() {return this->_updated;}

    // Set update state fo the particle
    void update(bool state) {this->_updated = state;}

    // Return private variable
    Type get(property p) {
        switch (p) {
            case POS :
                return this->_pos;
            case VEL :
                return this->_vel;
            case LOAD :
                return this->_load;
            default:
                break;
        }
    }

    float get_mass() {
        return this->_m;
    }

    // Set private variables from particle
    void set(property p, Type vec) {
        switch (p) {
            case POS :
                this->_pos = vec;
                break;
            case VEL :
                this->_vel = vec;
                break;
            case LOAD :
                this->_load = vec;
                break;
            default:
                break;
        }
    }

    void set(property p, float mass) {
        if (p == MASS) {
            this->_m = mass;
        }
    }

    void del_particle() {
        delete this->_pos;
        delete this->_vel;
        delete this->_load;
    }


private:
    // Set boolean to keep a state if updated already or not
    bool _updated = false;

    // Store position, velocity and load of the particle in space domain
    Type _pos;
    Type _vel;
    Type _load;

    // Store mass of the particle
    float _m = 0.0f;
};


#endif //PROJECT_PARTICLE_HPP
