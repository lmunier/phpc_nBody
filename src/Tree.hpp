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
        __host__ virtual ~AbstractType() {
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

        /**
         * Virtual function to be override in child class. Override by the child class.
         *
         * @return my_type which is the type of the class
         */
        __host__ __device__ virtual my_type get_type() = 0;

        /**
         * Return the mass of the given particle/cell.
         *
         * @return _m attribute of the given particle/cell
         */
        __host__ __device__ float get_mass() { return this->_m; }

        /**
         * Return the parent node of the given particle/cell.
         *
         * @return _parent attribute of the given particle/cell
         */
        __host__ __device__ AbstractType<Type>* get_parent() { return _parent; }

        /**
         * Virtual method to return the attribute vector array of pointer on the next particle/cells in the tree
         * data-structure. Override by the child class.
         *
         * @return _next attribute vector array of pointer on the next particle/cells in the tree data-structure
         */
        __host__ __device__ virtual thrust::device_vector<AbstractType < Type>* >* get_next() {}

        /**
         * Virtual method to get the value of an attribute of the given particle/cell. Override by the child class.
         *
         * @param p of the chosen attribute to be return
         * @return the value of the chosen attribute
         */
        __host__ __device__ virtual Type get(property p) = 0;

        /**
         * Set _m attribute value with the new mass
         *
         * @param mass new float value to set the attribute
         */
        __host__ __device__ void set_mass(float mass) { this->_m = mass; }

        /**
         * Set _parent attribute value with the new parent node
         *
         * @param parent new pointer on AbstractType<Type> value to set the attribute
         */
        __host__ __device__ void set_parent(AbstractType<Type>* parent) { this->_parent = parent; }

        /**
         * Virtual method to set attribute value of the chosen attribute for a given particle/cell. Override by the
         * child class.
         *
         * @param p enum value of the attribute
         * @param vec new vector value to set the attribute
         */
        __host__ __device__ virtual void set(property p, Type vec) {};

        /**
         * Find the index of the cell in the _next dynamic array where the particle should be stored.
         *
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param origin center of the parent cell
         * @param particle current particle to store
         * @return int index of the cell where the particle should be stored
         */
        virtual int find_cell_idx(Type origin) {};

        /**
         * Virtual method to compute the load vector between two elements to update the one passed in argument.
         * Override by the child class.
         *
         * @param particle where the load is applied
         */
        __device__ virtual void compute_load(AbstractType<Type>* particle) = 0;

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
        virtual bool is_out_boundaries() {};

        /**
         * Virtual method to update cell total mass and center of mass. Override by the child class.
         *
         * @param add boolean value to add or remove properties of the particle
         */
        virtual void update_cell(bool add) {};

        /**
         * Virtual method to update tree to change place of particles out of original boundaries. Override by the child
         * class.
         *
         * @return int value if the given particle/cell has to be delete (-1) or no (0)
         */
        virtual int update_tree() {};

        /**
         * Virtual method to subdivide a node in 2^NB_DIM sub-cells and fill _next attribute with a pointer to each
         * sub-cell. Override by the child class.
         */
        //__global__ virtual void subdivide_tree() {};

        /**
         * Virtual method store particle in the tree. Override by the child class.
         *
         * @param particle pointer on current particle to store
         * @param prev_part pointer on a particle to store after the current one
         */
        virtual void store_particle(AbstractType<Type>* particle, AbstractType<Type>* prev_part) {};

        /**
         * Virtual method to delete a given level of cells. Override by the child class.
         */
        virtual void del_level() {};

        /**
         * Virtual method to clear all elements in the attribute vector array of pointer on the next particle/cells in
         * the tree data-structure. Override by the child class.
         */
        virtual void clear_next() {};

    protected:
        bool _state = true;                         /**< state of the given particle/cell exist or not */
        float _m = 0.0f;                            /**< mass of the given particle/cell */
        AbstractType<Type>* _parent = nullptr;      /**< parent of the given particle/cell */
    };
}
#endif //PROJECT_TREE_HPP
