/**
 * @file Tree.hpp
 * @author Munier Louis
 * @date 27.05.19
 * @version 1.0
 *
 * Tree file to combine both Particle and Cell class in the same AbstractType. It is use to have vector array with both
 * elements.
 */

#ifndef PROJECT_TREE_HPP
#define PROJECT_TREE_HPP

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 * @include Vector.hpp custom library to have minimal vector implementation
 *
 * @namespace use std to simplify code implementation
 */
#include "constants.hpp"
#include "Vector.hpp"

using namespace std;

/**
 * Namespace of Tree to link classes.
 *
 * @namespace Tree
 */
namespace Tree {
    /**
     * @enum my_type to know which type of classe is in this namespace
     */
    enum my_type  {
        ParticleT,      /**< enum value for Particle type */
        CellT,          /**< enum value for Cell type */
    };

    /**
     * AbstractType class to link Particle class with Cell one.
     *
     * @class AbstractType
     */
    class AbstractType {
    public:
        /**
         * Virtual function to be override in child class.
         *
         * @return my_type which is the type of the class
         */
        virtual my_type get_type() = 0;

        /**
         * Destructor function to be override in child class. It serves as a model.
         */
        virtual ~AbstractType() = default;
    };

    /**
     * Class Particle, to store all the information of a given particle.
     *
     * @class Particle
     * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
     */
    template<typename Type>
    class Particle : public AbstractType {
    public:
        /**
         * Constructor to create a particle. All the member variables have defaults zero values if needed.
         *
         * @param mass of the given particle, in float
         * @param pos vector of the given particle
         */
        explicit Particle(float mass = 0.0f, Type pos = Type()) : _pos(pos), _m(mass) {}

        /**
         * Return the type of the given class.
         *
         * @return enum value ParticleT
         */
        my_type get_type() override { return ParticleT;}

        /**
         * Return boolean value if the given particle has already been updated.
         *
         * @return _updated boolean member variable
         */
        bool is_updated() { return this->_updated; }

        /**
         * Update the given value with the boolean state.
         *
         * @param state variable to fill boolean member variable _updated
         */
        void update(bool state) { this->_updated = state; }

        /**
         * Return one of the following properties of the given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the property
         * @return desired property of the given particle
         */
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

        /**
         * Return the mass of the given particle.
         *
         * @return _m member variable of the given particle
         */
        float get_mass() {
            return this->_m;
        }

        /**
         * Set member variable value of one of the following properties for a given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the property
         * @param vec new vector value to set the member variable
         */
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

        /**
         * Surcharge setting member variable value of one of the following properties for a given particle :
         * - mass
         *
         * @param p enum value of the property
         * @param mass new float value to set the member variable
         */
        void set(property p, float mass) {
            if (p == MASS) {
                this->_m = mass;
            }
        }

        /**
         * Delete a given particle.
         */
        void del_particle() {
            delete this->_pos;
            delete this->_vel;
            delete this->_load;
        }

        /**
         * Compute the load vector between two particle to update the one on the given particle.
         *
         * @param j_particle particle which influence the given particle
         */
        void compute_load(Particle<Type>* j_particle) {
            Type tmp = j_particle->get(POS) - this->get(POS);
            float norm = min(tmp.norm(), (float) EPSILON);

            this->_load = this->_load + tmp*(G*this->get_mass()*j_particle->get_mass()/tmp.norm());
        }

        /**
         * Update velocity and position of a given particle for a given load on it. Reset load after update.
         */
        void update_vel_pos() {
            Type new_velocity = this->get(LOAD) * (DELTA_T/this->get_mass()) + this->get(VEL);
            Type new_position = new_velocity * DELTA_T + this->get(POS);

            this->set(VEL, new_velocity);
            this->set(POS, new_position);
            this->set(LOAD, Type());
        }

    private:
        bool _updated = false;      /**< @var _updated, state of the given particle */

        Type _pos;                  /**< @var _pos, position of the given particle */
        Type _vel;                  /**< @var _vel, velocity of the given particle */
        Type _load;                 /**< @var _load, load on the given particle */

        float _m = 0.0f;            /**< @var _m, mass of the given particle */
    };

    /**
     * Class cell to store all the information of a given Cell.
     *
     * @class Cell
     * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
     */
    template<typename Type>
    class Cell : public AbstractType {
    public:
        /**
         * Constructor to create a Cell. All the member variables have defaults zero values if needed.
         *
         * @param center of the given cell
         * @param dim dimensions of the given cell
         * @param mass_pos center of mass of the given cell
         */
        explicit Cell(Type center = Type(), Type dim = Type(), Type mass_pos = Type()) {
            this->_center = center;
            this->_size = dim;
            this->_mass_pos = mass_pos;
        }

        /**
         * Return the type of the given class.
         *
         * @return enum value CellT
         */
        my_type get_type() override { return CellT;}

        /**
         * Subdivide a node in 2^NB_DIM sub-cells and fill _next member variable with a pointer to each sub-cell.
         */
        void subdivide_tree() {
            bool y = true, z = false;
            Type size = this->_size*0.5;

            /**< Loop to have 2^NB_DIM sub-cells */
            for(int n = 0; n < pow(2, NB_DIM); n++) {
                /**< Compute center of the cell */
                if (n%2 == 0)
                    y = !y;

                if (n == 4)
                    z = !z;

#if NB_DIM == DIM_2
                Type center = Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)));
#elif NB_DIM == DIM_3
                Type center = Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)), size.z*(0.5*pow(-1, z)));
#endif

                center = this->_center + center;

                /**< Fill _next vector array with each new sub-cell */
                auto next = new Cell<Type>(center, size, center);
                this->_next.push_back(next);
                next->_prev = this;
            }
        }

        /**
         * Delete a given particle.
         *
         * @param p pointer of the particle to be deleted
         */
        void del_element(Particle<Type>* p) {
            p->del_particle();

            for(auto it = this->_next.begin(); it != this->_next.end(); ++it) {
                if (it == p) {
                    (*it) = nullptr;
                    break;
                }
            }
        }

        float _m = 0.0f;                         /**< @var _m total mass of the particles in this cell */

        Type _center;                            /**< @var _center vector position of the center of the given cell */
        Type _size;                              /**< @var _size vector size of the given cell */
        Type _mass_pos;                          /**< @var _mass_pos vector position of the center of mass of the given
                                                  *                  cell */

        Cell* _prev;                             /**< @var _prev pointer on the parent cell*/
        vector< Tree::AbstractType*> _next{};    /**< @var _next vector list of the following cells or a particle */
    };
}
#endif //PROJECT_TREE_HPP
