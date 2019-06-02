/**
 * @file Particle.hpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
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
#include "Vector.hpp"
#include "Tree.hpp"

/**
 * @namespace std to simplify code implementation
 */
using namespace std;

/**
 * Namespace of Tree to link classes.
 *
 * @namespace Tree
 */
namespace Tree {
    /**
     * Class Particle, to store all the information of a given particle.
     *
     * @class Particle
     * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
     */
    template<typename Type>
    class Particle : public AbstractType<Type> {
    public:
        /**
         * Constructor to create a particle. All the member variables have defaults zero values if needed.
         *
         * @param mass of the given particle, in float
         * @param pos vector of the given particle
         */
        explicit Particle(float mass = 0.0f, Type pos = Type()) : _pos(pos) {
            this->set_mass(mass);
        }

        /**
         * Destructor to safely delete a given particle.
         */
        ~Particle() = default;

        /**
         * Return one of the following attribute of the given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the attribute
         * @return desired attribute of the given particle
         */
        Type get(property p) override {
            switch (p) {
                case POS :
                    return this->_pos;
                case VEL :
                    return this->_vel;
                case LOAD :
                    return this->_load;
                default:
                    return Type();
                    break;
            }
        }

        /**
         * Set attribute value of one of the following attribute for a given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the attribute
         * @param vec new vector value to set the attribute
         */
        void set(property p, Type vec) override {
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

        /**
         * Update velocity and position of a given particle for a given load on it. Reset load after update.
         */
        void update_vel_pos() override {
            Type new_velocity = this->get(LOAD) * (DELTA_T / this->get_mass()) + this->get(VEL);
            Type new_position = new_velocity * DELTA_T + this->get(POS);

            this->set(VEL, new_velocity);
            this->set(POS, new_position);
            this->set(LOAD, Type());
        }

        /**
         * Test if a given particle is out of the overall area and return true if it is the case, false otherwise.
         *
         * @return true if particle out of its boundaries, false otherwise
         */
        bool is_out_boundaries(Type dim) {
            Type pos = this->get(POS);

            if (abs(pos.x) > 0.5 * dim.x)
                return true;

            if (abs(pos.y) > 0.5 * dim.y)
                return true;

#if NB_DIM == DIM_3
            if (abs(pos.z) > 0.5 * dim.z)
                return true;
#endif
            return false;
        }

    private:
        Type _pos;                  /**< position of the given particle */
        Type _vel;                  /**< velocity of the given particle */
        Type _load;                 /**< load on the given particle */
    };
}


#endif //PHPC_NBODY_PARTICLE_HPP
