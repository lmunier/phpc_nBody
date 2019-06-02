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
    template<typename Type>
    class AbstractType {
    public:
        /**
         * Destructor function to safely delete all pointers in this class and set state of AbstractType to false.
         */
        __host__ __device__ virtual ~AbstractType() {
            this->_state = false;
            delete this->_parent;
        };

        __host__ void *operator new(size_t len) {
            void *ptr;
            cudaMallocManaged(&ptr, len);
            cudaDeviceSynchronize();
            return ptr;
        }

        __host__ void operator delete(void *ptr) {
            cudaDeviceSynchronize();
            cudaFree(ptr);
        }

        /**
         * State of the cell, if it is already deleted or not.
         *
         * @return state of the particle/cell already deleted (false) or not (true)
         */
        __host__ __device__ bool get_state() { return this->_state; }

    protected:
        bool _state = true;                         /**< state of the given particle/cell exist or not */
        float _m = 0.0f;                            /**< mass of the given particle/cell */
        AbstractType<Type>* _parent = nullptr;      /**< parent of the given particle/cell */
    };
}
#endif //PROJECT_TREE_HPP
