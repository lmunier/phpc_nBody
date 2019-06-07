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
         * Return the type of the given class.
         *
         * @return enum value ParticleT
         */
        my_type get_type() override { return ParticleT; }

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
         * Find the index of the cell in the _next dynamic array where the particle should be stored.
         *
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param origin center of the parent cell
         * @param particle current particle to store
         * @return int index of the cell where the particle should be stored
         */
        int find_cell_idx(Type origin) override {
            int idx = 0;
            Type tmp_vec = this->get(POS) - origin;

            if(tmp_vec.x > 0)
                idx += 1;

            if(tmp_vec.y < 0)
                idx += 2;

#if NB_DIM == DIM_3
            if(tmp_vec.z < 0)
                idx += 4;
#endif

            return idx;
        }

        /**
         * Compute the load vector between two particles to update the one passed in argument. The load is
         * the General gravity one.
         *
         * @param particle where the load is applied
         */
        void compute_load(AbstractType<Type> *particle) override {
            Type tmp = this->get(POS) - particle->get(POS);
            float d = max(tmp.norm(), EPSILON);
            particle->set(LOAD, particle->get(LOAD) + tmp * (G * particle->get_mass() * this->get_mass()) / d);
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
         * Test if a given particle is out of its cell and return true if it is the case, false otherwise.
         *
         * @return true if particle out of its boundaries, false otherwise
         */
        bool is_out_boundaries() {
            auto parent = this->get_parent();

            if (parent == nullptr)
                return true;

            Type distance = this->get(POS) - parent->get(CENTER);
            Type cell_size = parent->get(DIM);

            if (2 * abs(distance.x) > cell_size.x)
                return true;

            if (2 * abs(distance.y) > cell_size.y)
                return true;

#if NB_DIM == DIM_3
            if (2 * abs(distance.z) > cell_size.z)
                return true;
#endif
            return false;
        }

        /**
         * Update cell total mass and center of mass.
         *
         * @param add boolean value to add (true) or remove (false) attribute values of the particle
         */
        void update_cell(bool add) override {
            auto head = this->get_parent();

            float mass_this = this->get_mass();
            float mass_head = head->get_mass();

            if (add) {
                float mass_tot = mass_head + mass_this;

                head->set(MASS_POS, (head->get(MASS_POS) * mass_head + this->get(POS) * mass_this) / mass_tot);
                head->set_mass(mass_tot);
            } else {
                float mass_tot = max(mass_head - mass_this, (float) EPSILON);

                /** Avoid having a division by zero */
                if (mass_tot == 0.0f)
                    head->set(MASS_POS, head->get(CENTER));
                else
                    head->set(MASS_POS, (head->get(MASS_POS) * mass_head - this->get(POS) * mass_this) / mass_tot);

                head->set_mass(mass_tot);
            }
        }

        /**
         * Update tree to change place of particles out of original boundaries. Override by the child
         * class.
         *
         * @return int value if the given particle/cell has to be delete (-1) or no (0)
         */
        int update_tree() override {
            bool first = true;

            /** continue to go up in level to find right cell for our particle */
            while (this->is_out_boundaries()) {
                int nb_particles = 0;
                auto parent = this->get_parent();

                if (parent == nullptr)
                    return -1;

                /** delete pointer from cell to particle */
                if (first) {
                    parent->clear_next();
                    first = false;
                } else {
                    /** Check for empty level */
                    for (auto n : parent->get_next()) {
                        if (!(n->get_next().empty())) {
                            nb_particles++;
                            break;
                        }
                    }

                    /** delete empty level of the parent node */
                    if (nb_particles == 0) {
                        parent->del_level();
                        return -2;
                    }
                }

                this->update_cell(false);
                this->set_parent(parent->get_parent());
            }

            this->update_cell(false);
            this->get_parent()->store_particle(this, nullptr);
            return 0;
        }

    private:
        Type _pos;                  /**< position of the given particle */
        Type _vel;                  /**< velocity of the given particle */
        Type _load;                 /**< load on the given particle */
    };
}


#endif //PHPC_NBODY_PARTICLE_HPP
