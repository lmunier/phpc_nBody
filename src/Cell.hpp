/**
 * @file Particle.hpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
 *
 * Cell file to implement a class named cell,
 */

#ifndef PHPC_NBODY_CELL_HPP
#define PHPC_NBODY_CELL_HPP

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
     * Class cell to store all the information of a given Cell.
     *
     * @class Cell
     * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
     */
    template<typename Type>
    class Cell : public AbstractType<Type> {
    public:
        /**
         * Constructor to create a Cell. All the member variables have defaults zero values if needed.
         *
         * @param center center of the given cell
         * @param dim dimensions of the given cell
         * @param mass_pos center of mass of the given cell
         */
        __host__ explicit Cell(Type center = Type(), Type dim = Type(), Type mass_pos = Type()) :
                _center(center), _size(dim), _mass_pos(mass_pos) {}

        /**
         * Destructor to safely delete all pointer in the Cell element.
         */
        __host__ __device__ ~Cell() override {
            this->set_parent(nullptr);
            delete this->get_parent();

            for (auto n : this->_next)
                delete n;

            this->_next.clear();
        }

        /**
         * Return the type of the given class.
         *
         * @return enum value CellT
         */
        __host__ __device__ my_type get_type() override { return CellT; }

        /**
         * Return the attribute vector array of pointer on the next particle/cells in the tree data-structure. Override
         * by the child class.
         *
         * @return _next attribute vector array of pointer on the next particle/cells in the tree data-structure
         */
        __device__ vector<AbstractType < Type>* > get_next() override { return _next; }

        /**
         * Return one of the following attribute of the given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the attribute
         * @return desired attribute of the given particle
         */
        __host__ __device__ Type get(property p) override {
            switch (p) {
                case CENTER :
                    return this->_center;
                case DIM :
                    return this->_size;
                case MASS_POS :
                    return this->_mass_pos;
                default:
                    break;
            }
        }

        /**
         * Set attribute value of one of the following properties for a given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the attribute
         * @param vec new vector value to set the attribute
         */
        __host__ __device__ void set(property p, Type vec) override {
            switch (p) {
                case CENTER :
                    this->_center = vec;
                    break;
                case DIM :
                    this->_size = vec;
                    break;
                case MASS_POS :
                    this->_mass_pos = vec;
                    break;
                default:
                    break;
            }
        }

        /**
         * Store particle in the tree.
         *
         * @param particle pointer on current particle to store
         * @param prev_part pointer on a particle to store after the current one
         */
        // void store_particle(AbstractType <Type> *particle, AbstractType <Type> *prev_part) override {
        //     /** If there is no next element */
        //     if (this->_next.empty()) {
        //         this->_next.push_back(particle);

        //         particle->set_parent(this);
        //         particle->update_cell(true);

        //         if (prev_part != nullptr) {
        //             particle = prev_part;
        //             prev_part = nullptr;

        //             this->get_parent()->store_particle(particle, prev_part);
        //         }
        //     } /** If next element is a particle */
        //     else if (this->_next[0]->get_type() == ParticleT) {
        //         prev_part = this->_next[0];

        //         /** clear precedent pointers */
        //         this->_next.clear();
        //         this->set_mass(0.0f);
        //         this->set(MASS_POS, Type());
        //         prev_part->set_parent(nullptr);

        //         this->subdivide_tree();
        //         this->store_particle(particle, prev_part);
        //     } /** If next element is a cell */
        //     else {
        //         int cell_idx_1 = particle->find_cell_idx(this->_center);

        //         particle->set_parent(this);
        //         particle->update_cell(true);
        //         particle->set_parent(nullptr);

        //         this->get_next()[cell_idx_1]->store_particle(particle, prev_part);
        //     }
        // }

        /**
         * Compute the load vector between two particles to update the one passed in argument. Two different types of
         * load are available by changing a define in the constants.hpp file :
         * - General gravity load
         * - Lennard Jones potential
         *
         * @param particle where the load is applied
         */
        __device__ void compute_load(AbstractType <Type> *particle) override {
#if LOAD_TYPE == 0
            Type tmp = this->get(MASS_POS) - particle->get(POS);
            float d = max(tmp.norm(), EPSILON);
            particle->set(LOAD, particle->get(LOAD) + tmp * (G * particle->get_mass() * this->get_mass()) / d);
#elif LOAD_TYPE == 1
            Type tmp = particle->get(POS) - this->get(MASS_POS);
            float d = tmp.norm();
            particle->set(LOAD, particle->get(LOAD) + tmp * (48 / pow(d, 2)) * (1 / pow(d, 12) - 0.5f * pow(d, 6)));
#endif
        }

//         /**
//          * Delete a given level of cells
//          */
//         void del_level() override {
//             while (!this->_next.empty()) {
//                 auto next = this->_next.back();

//                 next->set_parent(nullptr);
//                 delete next->get_parent();

//                 this->_next.pop_back();
//                 delete next;
//             }

//             this->_next.clear();
//         }

//         /**
//          * Clear all elements in the attribute vector array of pointer on the next particle/cells in the tree
//          * data-structure. Override by the child class.
//          */
//         void clear_next() override {
//             this->_next.clear();
//         }

    private:
        Type _center;                            /**< vector position of the center of the given cell */
        Type _size;                              /**< vector size of the given cell */
        Type _mass_pos;                          /**< vector position of the center of mass of the given cell */

        vector<AbstractType < Type>*> _next{};    /**< vector list of the following cells or a particle */
    };
}

#endif //PHPC_NBODY_CELL_HPP
