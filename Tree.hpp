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
    template<typename Type>
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

        /**
         * Return the mass of the given particle/cell.
         *
         * @return _m member variable of the given particle
         */
        float get_mass() {
            return this->_m;
        }

        /**
         * Surcharge setting member variable value of one of the following properties for a given particle :
         * - mass
         *
         * @param p enum value of the property
         * @param mass new float value to set the member variable
         */
        void set_mass(float mass) {
            this->_m = mass;
        }

        virtual Type get(property p) = 0;
        virtual AbstractType<Type>* get_parent() {};
        virtual void set(property p, Type vec) {};
        virtual void set_parent(AbstractType<Type>*) {};
        virtual void compute_load(AbstractType<Type>* particle) = 0;
        virtual void update_vel_pos() {};
        virtual bool is_out_boundaries() {};
        virtual void update_cell(bool add) {};
        virtual int update_tree() {};

        virtual void subdivide_tree() {};
        virtual void store_particle(AbstractType<Type>* particle, vector< AbstractType<Type>* > *list_p) {};
        virtual void del_level() {};
        virtual AbstractType<Type>* get_prev() {};

        virtual vector< AbstractType<Type>* > get_next() {};

        float _m = 0.0f;            /**< @var mass of the given particle/cell */
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
         * Destructor to correctly delete a given particle.
         */
        ~Particle() override {
            delete this->_parent;
        }
        /**
         * Return the type of the given class.
         *
         * @return enum value ParticleT
         */
        my_type get_type() override { return ParticleT;}

        /**
         * Return one of the following properties of the given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the property
         * @return desired property of the given particle
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

        AbstractType<Type>* get_parent() {
            return this->_parent;
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

        void set_parent(AbstractType<Type>* parent) {
            this->_parent = parent;
        }

        /**
         * Compute the load vector between two particles to update the one passed in argument.
         *
         * @param particle where the load is applied
         */
        void compute_load(AbstractType<Type>* particle) override {
            Type tmp = this->get(POS) - particle->get(POS);
            float norm = min(tmp.norm(), (float) EPSILON);

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
         * Test if a particle is out of its cell and return true if it is the case, false otherwise.
         *
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param particle pointer on the tested particle
         * @return true if particle out of its boundaries, false otherwise
         */
        bool is_out_boundaries() {
            Type position = this->get(POS);
            Type center = this->_parent->get(CENTER);
            Type cell_size = this->_parent->get(DIM);

            Type distance = position - center;

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
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param particle that need to update its cell
         * @param add boolean value to add or remove properties of the particle
         */
        void update_cell(bool add){
            auto head = this->_parent;

            if (add) {
                head->set(MASS_POS, (head->get(MASS_POS) * head->get_mass() + this->get(POS) * this->get_mass())
                                    / (head->get_mass() + this->get_mass()));
                head->set_mass(head->get_mass() + this->get_mass());
            } else {
                head->set(MASS_POS, (head->get(MASS_POS) * head->get_mass() - this->get(POS) * this->get_mass())
                                    / (head->get_mass() - this->get_mass()));
                head->set_mass(head->get_mass() - this->get_mass());
            }
        }

        /**
         * Update tree to change place of particles out of original boundaries.
         *
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param particle that change its place
         */
        int update_tree() override {
            int nb_particles = 0;
            auto parent = this->_parent;

            /**< delete pointer from cell to particle */
            parent->get_next().clear();

            /**< parent of the particle go back from one level */
            this->_parent = parent->get_prev();
            parent = parent->get_prev();

            if (this->is_out_boundaries()) {
                for (auto it = parent->get_next().begin(); it != parent->get_next().end(); ++it) {
                    if (!(*it)->get_next().empty())
                        nb_particles++;
                }

                /**< delete empty level of the parent node */
                if (nb_particles == 0)
                    this->_parent->del_level();
            }

            /**< continue to go up in level to find right cell for our particle */
            while(this->is_out_boundaries()) {
                this->update_cell(false);
                this->_parent = this->_parent->get_prev();

                if (this->_parent == nullptr) {
                    delete this;
                    return -1;
                }
            }

            return 0;
        }


    private:
        AbstractType<Type>* _parent = nullptr;

        Type _pos;                  /**< @var position of the given particle */
        Type _vel;                  /**< @var velocity of the given particle */
        Type _load;                 /**< @var load on the given particle */
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
         * @param center of the given cell
         * @param dim dimensions of the given cell
         * @param mass_pos center of mass of the given cell
         */
        explicit Cell(Type center = Type(), Type dim = Type(), Type mass_pos = Type()) :
            _center(center), _size(dim), _mass_pos(mass_pos) {}

        /**
         * Destructor to correctly delete all pointer in the Cell element.
         */
        ~Cell() override{
            delete this->_prev;

            while(!this->_next.empty()) {
                delete this->_next.back();
                this->_next.pop_back();
            }
        }

        /**
         * Return the type of the given class.
         *
         * @return enum value CellT
         */
        my_type get_type() override { return CellT;}

        /**
         * Return one of the following properties of the given particle :
         * - position
         * - velocity
         * - load
         *
         * @param p enum value of the property
         * @return desired property of the given particle
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

        AbstractType<Type>* get_prev() override {
            return this->_prev;
        }

        vector< AbstractType<Type>* > get_next() override {
            return this->_next;
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
         * Subdivide a node in 2^NB_DIM sub-cells and fill _next member variable with a pointer to each sub-cell.
         */
        void subdivide_tree() override {
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
         * Store particle in the tree.
         *
         * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
         * @param head pointer on current cell to store the particle
         * @param particle pointer on current particle to store
         * @param list_p pointer on a list of pointers on all the particles to store after the current one
         */
        void store_particle(AbstractType<Type>* particle, vector< AbstractType<Type>* > *list_p) override {
            int cell_idx_1 = find_cell_idx(this->_center, particle->get(POS));

            if (list_p == nullptr) {
                vector< AbstractType<Type>* > other_particle;
                list_p = &other_particle;
            }

            if (this->_next.empty()) {
                this->_next.push_back(particle);
                particle->set_parent(this);
                particle->update_cell(true);

                Cell<Type>* previous;
                do {
                    previous = this->_prev;
                    if (previous == nullptr)
                        break;

                    particle->update_cell(true);
                } while (previous->_prev != nullptr);

                if (!list_p->empty()) {
                    particle = list_p->back();
                    list_p->pop_back();

                    this->store_particle(particle, list_p);
                }
            } else if (this->_next[0]->get_type() == ParticleT) {
                list_p->push_back((Particle<Type>*) this->_next[0]);
                this->_next[0]->update_cell(false);

                /**< clear precedent pointers */
                this->_next.clear();
                particle->set_parent(nullptr);

                this->subdivide_tree();
                this->store_particle(particle, list_p);
            } else {
                this->get_next()[cell_idx_1]->store_particle(particle, list_p);
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
         * Delete a given particle.
         *
         * @param p pointer of the particle to be deleted
         */
        void del_level() override {
            while(!this->_next.empty()) {
                delete this->_next.back();
                this->_next.pop_back();
            }
        }

        Type _center;                            /**< @var vector position of the center of the given cell */
        Type _size;                              /**< @var vector size of the given cell */
        Type _mass_pos;                          /**< @var vector position of the center of mass of the given cell */

        Cell* _prev;                             /**< @var pointer on the parent cell*/
        vector< AbstractType<Type>*> _next{};    /**< @var vector list of the following cells or a particle */
    };
}
#endif //PROJECT_TREE_HPP
