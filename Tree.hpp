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
         * Virtual function to be override in child class. Override by the child class.
         *
         * @return my_type which is the type of the class
         */
        virtual my_type get_type() = 0;

        /**
         * Destructor function to safely delete all pointers in this class and set state of AbstractType to false.
         */
        virtual ~AbstractType() {
            this->_state = false;
            delete this->_parent;
        };

        /**
         * Return the mass of the given particle/cell.
         *
         * @return _m attribute of the given particle/cell
         */
        float get_mass() { return this->_m; }

        /**
         * Return the parent node of the given particle/cell.
         *
         * @return _parent attribute of the given particle/cell
         */
        AbstractType<Type>* get_parent() { return this->_parent; }

        /**
         * Virtual method to return the initialization state of the given particle/cell. Override by the child class.
         *
         * @return _init attribute of the given particle/cell
         */
        virtual bool get_init() {};

        /**
         * Set _m attribute value with the new mass
         *
         * @param mass new float value to set the attribute
         */
        void set_mass(float mass) { this->_m = mass; }

        /**
         * Set _parent attribute value with the new parent node
         *
         * @param parent new pointer on AbstractType<Type> value to set the attribute
         */
        void set_parent(AbstractType<Type>* parent) { this->_parent = parent; }

        /**
         * Virtual method to set the new value of _init attribute of the given particle/cell. Override by the child
         * class.
         *
         * @param init value to set _init attribute
         */
        virtual void set_init(bool init) {};

        /**
         * Virtual method to get the value of an attribute of the given particle/cell. Override by the child class.
         *
         * @param p of the chosen attribute to be return
         * @return the value of the chosen attribute
         */
        virtual Type get(property p) = 0;

        /**
         * Virtual method to set attribute value of the chosen attribute for a given particle/cell. Override by the
         * child class.
         *
         * @param p enum value of the attribute
         * @param vec new vector value to set the attribute
         */
        virtual void set(property p, Type vec) {};

        /**
         * Virtual method to compute the load vector between two elements to update the one passed in argument.
         * Override by the child class.
         *
         * @param particle where the load is applied
         */
        virtual void compute_load(AbstractType<Type>* particle) = 0;

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
        virtual void subdivide_tree() {};

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
         * Virtual method to return the attribute vector array of pointer on the next particle/cells in the tree
         * data-structure. Override by the child class.
         *
         * @return _next attribute vector array of pointer on the next particle/cells in the tree data-structure
         */
        virtual vector< AbstractType<Type>* > get_next() {};

        /**
         * Virtual method to clear all elements in the attribute vector array of pointer on the next particle/cells in
         * the tree data-structure. Override by the child class.
         */
        virtual void clear_next() {};

        /**
         * Virtual method to set the _prev AbstracteType pointer attribute to the previous particle/cell with the new
         * previous value. Override by the child class.
         *
         * @param previous pointer attribute to the previous particle/cell with the new previous value
         */
        virtual void set_prev(AbstractType<Type>* previous) {};

        bool _state = true;                         /**< state of the given particle/cell exist or not */
        float _m = 0.0f;                            /**< mass of the given particle/cell */
        AbstractType<Type>* _parent = nullptr;      /**< parent of the given particle/cell */
    };

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
        my_type get_type() override { return ParticleT;}

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
         * Set the new value of _init attribute of the given particle/cell. Override by the child class.
         *
         * @param init value to set _init attribute
         */
        bool get_init() override { return this->_init; }

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
         * Set the new value of _init attribute of the given particle/cell. Override by the child class.
         *
         * @param init value to set _init attribute
         */
        void set_init(bool init) override {this->_init = init; }

        /**
         * Compute the load vector between two particles to update the one passed in argument.
         *
         * @param particle where the load is applied
         */
        void compute_load(AbstractType<Type>* particle) override {
            Type tmp = this->get(POS) - particle->get(POS);
            float norm = max(tmp.norm(), (float) EPSILON);

            particle->set(LOAD, particle->get(LOAD) + tmp*(G*particle->get_mass()*this->get_mass())/tmp.norm());
        }

        /**
         * Update velocity and position of a given particle for a given load on it. Reset load after update.
         */
        void update_vel_pos() override {
            Type new_velocity = this->get(LOAD) * (DELTA_T/this->get_mass()) + this->get(VEL);
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
            if (this->get_parent() == nullptr)
                return true;

            Type distance = this->get(POS) - this->get_parent()->get(CENTER);
            Type cell_size = this->get_parent()->get(DIM);

            if (2*abs(distance.x) > cell_size.x)
                return true;

            if (2*abs(distance.y) > cell_size.y)
                return true;

#if NB_DIM == DIM_3
            if (2*abs(distance.z) > cell_size.z)
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

            if (add) {
                float mass_tot = head->get_mass() + this->get_mass();

                head->set(MASS_POS, (head->get(MASS_POS) * head->get_mass() + this->get(POS) * this->get_mass())
                                    / mass_tot);
                head->set_mass(mass_tot);
            } else {
                float mass_tot = max(head->get_mass() - this->get_mass(), (float) EPSILON);

                /** Avoid having a division by zero */
                if (mass_tot == 0.0f)
                    head->set(MASS_POS, head->get(CENTER));
                else
                    head->set(MASS_POS, (head->get(MASS_POS) * head->get_mass() - this->get(POS) * this->get_mass())
                                        / mass_tot);

                head->set_mass(mass_tot);
            }

            if (this->_init)
                return;
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

                if (this->get_parent() == nullptr)
                    return -1;

                /** delete pointer from cell to particle */
                if (first) {
                    this->get_parent()->clear_next();
                    first = false;
                } else {
                    /** Check for empty level */
                    for (auto n : this->get_parent()->get_next()) {
                        if (!(n->get_next().empty())) {
                            nb_particles++;
                            break;
                        }
                    }

                    /** delete empty level of the parent node */
                    if (nb_particles == 0)
                        this->get_parent()->del_level();
                }

                this->update_cell(false);
                this->set_parent(this->get_parent()->get_parent());
            }

            this->update_cell(false);
            this->get_parent()->store_particle(this, nullptr);
            return 0;
        }

    private:
        bool _init = false;         /**< initialization state of the given particle */

        Type _pos;                  /**< position of the given particle */
        Type _vel;                  /**< velocity of the given particle */
        Type _load;                 /**< load on the given particle */
    };

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
        explicit Cell(Type center = Type(), Type dim = Type(), Type mass_pos = Type()) :
                _center(center), _size(dim), _mass_pos(mass_pos) {}

        /**
         * Destructor to safely delete all pointer in the Cell element.
         */
        ~Cell() {
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
        my_type get_type() override { return CellT;}

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
         * Return the attribute vector array of pointer on the next particle/cells in the tree data-structure. Override
         * by the child class.
         *
         * @return _next attribute vector array of pointer on the next particle/cells in the tree data-structure
         */
        vector< AbstractType<Type>* > get_next() override {
            return this->_next;
        }

        /**
         * Clear all elements in the attribute vector array of pointer on the next particle/cells in the tree
         * data-structure. Override by the child class.
         */
        void clear_next() override {
            this->_next.clear();
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
        void set(property p, Type vec) override {
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
         * Subdivide a node in 2^NB_DIM sub-cells and fill _next attribute vector array with a pointer to each sub-cell.
         */
        void subdivide_tree() override {
            bool y = true, z = false;
            Type size = this->_size*0.5;

            /** Loop to have 2^NB_DIM sub-cells */
            for(int n = 0; n < pow(2, NB_DIM); n++) {
                /** Compute center of the cell */
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

                /** Fill _next vector array with each new sub-cell */
                auto next = new Cell<Type>(center, size, center);
                this->_next.push_back(next);
                next->set_parent(this);
            }
        }

        /**
         * Store particle in the tree.
         *
         * @param particle pointer on current particle to store
         * @param prev_part pointer on a particle to store after the current one
         */
        void store_particle(AbstractType<Type>* particle, AbstractType<Type>* prev_part) override {
            /** If there is no next element */
            if (this->_next.empty()) {
                this->_next.push_back(particle);

                particle->set_parent(this);
                particle->update_cell(true);

                if (prev_part != nullptr) {
                    particle = prev_part;
                    prev_part = nullptr;

                    this->get_parent()->store_particle(particle, prev_part);
                }
            } /** If next element is a particle */
            else if (this->_next[0]->get_type() == ParticleT) {
                prev_part = this->_next[0];

                /** clear precedent pointers */
                this->_next.clear();
                this->set_mass(0.0f);
                this->set(MASS_POS, Type());
                prev_part->set_parent(nullptr);

                this->subdivide_tree();
                this->store_particle(particle, prev_part);
            } /** If next element is a cell */
            else {
                int cell_idx_1 = find_cell_idx(this->_center, particle->get(POS));

                particle->set_parent(this);
                particle->update_cell(true);
                particle->set_parent(nullptr);

                this->get_next()[cell_idx_1]->store_particle(particle, prev_part);
            }
        }

        /**
         * Compute the load vector between a particle and a particle cluster represented by a center of mass.
         *
         * @param particle where the load is applied
         */
        void compute_load(AbstractType<Type>* particle) override {
            Type tmp = this->_mass_pos - particle->get(POS);
            float norm = tmp.norm();

            particle->set(LOAD, particle->get(LOAD) + tmp*(G*particle->get_mass()*this->_m)/tmp.norm());
        }

        /**
         * Delete a given level of cells
         */
        void del_level() override {
            while (!this->_next.empty()) {
                auto next = this->_next.back();

                next->set_parent(nullptr);
                delete next->get_parent();

                this->_next.pop_back();
                delete next;
            }

            this->_next.clear();
        }

        Type _center;                            /**< vector position of the center of the given cell */
        Type _size;                              /**< vector size of the given cell */
        Type _mass_pos;                          /**< vector position of the center of mass of the given cell */

        vector< AbstractType<Type>*> _next{};    /**< vector list of the following cells or a particle */
    };
}
#endif //PROJECT_TREE_HPP
