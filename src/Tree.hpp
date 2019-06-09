/**
 * @file Tree.hpp
 * @author Munier Louis
 * @date 30.05.19
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
     * AbstractType class to link Particle class with Cell one.
     *
     * @class AbstractType
     */
    template<typename Type>
    class AbstractType {
    public:
        /**
         * Destructor function to safely delete all pointers in this class and set state of AbstractType to false.
         */
        virtual ~AbstractType() = default;

        /**
         * Return the mass of the given particle/cell.
         *
         * @return _m attribute of the given particle/cell
         */
        float get_mass() { return this->_m; }

        /**
         * Virtual method to get the value of an attribute of the given particle/cell. Override by the child class.
         *
         * @param p of the chosen attribute to be return
         * @return the value of the chosen attribute
         */
        virtual Type get(property p) = 0;

        /**
         * Set _m attribute value with the new mass
         *
         * @param mass new float value to set the attribute
         */
        void set_mass(float mass) { this->_m = mass; }

        /**
         * Virtual method to set attribute value of the chosen attribute for a given particle/cell. Override by the
         * child class.
         *
         * @param p enum value of the attribute
         * @param vec new vector value to set the attribute
         */
        virtual void set(property p, Type vec) {};

        /**
         * Virtual method to update velocity and position of a given particle for a given load on it. Reset load after
         * update. Override by the child
         */
        virtual void update_vel_pos() {};

        /**
         * Virtual method to test if a particle is out of its cell and return true if it is the case, false otherwise.
         * Override by the child class.
         *
         * @return true if particle out of its boundaries, false otherwise
         */
        virtual bool is_out_boundaries(Type dim) { return false; };

    protected:
        float _m = 0.0f;                            /**< mass of the given particle/cell */
    };
}
#endif //PROJECT_TREE_HPP
