/**
 * @file Tree.cpp
 * @author Munier Louis
 * @date 27.05.19
 * @version 1.0
 *
 * Tree file to combine both Particle and Cell class in the same AbstractType. It is use to have vector array with both
 * elements. It also implement some useful functions related to tree data-structure.
 */

#include "Tree.hpp"

/**
 * Test if a particle is out of its cell and return true if it is the case, false otherwise.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle pointer on the tested particle
 * @return true if particle out of its boundaries, false otherwise
 */
template<typename Type>
bool Tree::Particle<Type>::is_out_boundaries() {
    Type distance = this->get(POS) - dynamic_cast<Cell<Type> *>(this->parent)->_center;
    Type cell_size = this->parent->_size;

    if (2*abs(distance.x) < cell_size.x)
        return true;
    else if (2*abs(distance.y) < cell_size.y)
        return true;
#if NB_DIM == DIM_3
    if (2*abs(distance.z) < cell_size.z)
        return true;
#endif

    return false;
}

/**
 * Generate NB_PARTICLES particles and call function to store them.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer to initial Cell of the tree to construct
 * @param vec vector to give the Type of the vector to the template function
 * @param list_particles pointer on a list of all the pointers on generated particles (if PRINT defined) to generate csv
 */
template <typename Type>
void Tree::Cell<Type>::generate_data(vector< Particle<Type>* > &list_particles) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5*this->_size.x;
    float y_rnd, p_y = 0.5*this->_size.y;

    vector< Particle<Type>* > other_particle;

    random_device rd;
    uniform_real_distribution<float> dist_m(0, p_m);
    uniform_real_distribution<float> dist_x(-p_x, p_x);
    uniform_real_distribution<float> dist_y(-p_y, p_y);

    for (int i = 0; i < NB_PARTICLES; i++) {
        /**< Generate random coordinate for a new particle and shift if it stays in boundaries to stress application */
        x_rnd = dist_x(rd);
        if(fabs(x_rnd - SHIFT) < p_x)
            x_rnd -= SHIFT;

        y_rnd = dist_y(rd);
        if(fabs(y_rnd - SHIFT) < p_y)
            y_rnd -= SHIFT;

#if NB_DIM == DIM_3
        float z_rnd, p_z = 0.5*this->_size.z;
        uniform_real_distribution<float> dist_z(-p_z, p_z);

        z_rnd = dist_z(rd);
        if(fabs(z_rnd - SHIFT) < p_z)
            z_rnd -= SHIFT;
#endif

        /**< Initialize quadtree/octree */
        if(this->_prev == this) {
            this->subdivide_tree();
            this->_prev = nullptr;
        }

        /**< Create new particle */
#if NB_DIM == DIM_2
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

        list_particles.push_back(new_particle);
        this->store_particle(new_particle, &other_particle);
    }
}

/**
 * Subdivide a node in 2^NB_DIM sub-cells and fill _next member variable with a pointer to each sub-cell.
 */
template <typename Type>
void Tree::Cell<Type>::subdivide_tree() {
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
template <typename Type>
void Tree::Cell<Type>::store_particle(Particle<Type>* particle, vector< Particle<Type>* >* list_p) {
    int cell_idx_1 = find_cell_idx(this->_center, particle->get(POS));

    if (this->_next.empty()) {
        this->_next.push_back(particle);
        this->_mass_pos = (this->_mass_pos * this->_m + particle->get(POS) * particle->get_mass()) *
                          ( 1 / (this->_m + particle->get_mass()));
        this->_m += particle->get_mass();

        particle->parent = this;

        do {
            this = this->_prev;
            this->_mass_pos = (this->_mass_pos * this->_m + particle->get(POS) * particle->get_mass()) *
                              ( 1 / (this->_m + particle->get_mass()));
            this->_m += particle->get_mass();

        } while (this->_prev != nullptr);

        if (!list_p->empty()) {
            particle = list_p->back();
            list_p->pop_back();

            store_particle(this, particle, list_p);
        }
    } else if (this->_next[0]->get_type() == ParticleT) {
        list_p->push_back((Particle<Type>*) this->_next[0]);
        this->_mass_pos = (this->_mass_pos * this->_m - particle->get(POS) * particle->get_mass()) *
                          ( 1 / (this->_m - particle->get_mass()));
        this->_m -= particle->get_mass();

        /**< clear precedent pointers */
        this->_next.clear();
        particle->parent = nullptr;

        this->subdivide_tree();
        store_particle(this, particle, list_p);
    } else {
        store_particle((Cell<Type>*) this->_next[cell_idx_1], particle, list_p);
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
template <typename Type>
int Tree::Cell<Type>::find_cell_idx(Type origin, Type particle) {
    int idx = 0;
    Type tmp_vec = particle - origin;

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
 * Update the load applied to a particle.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param head pointer to the current cell of the tree
 * @param part_loaded pointer to the current particle for which the load is computed
 */
/*template <typename Type>
void Tree::update_load(Cell<Type> *head, Particle<Type> *part_loaded) {
    for (auto it = head->_next.begin(); it != head->_next.end(); ++it) {
        if ((*it) == nullptr)
            return;
        else if ((*it)->get_type() == ParticleT) {
            if ((*it) != part_loaded)
                dynamic_cast<Particle<Type> *>(*it)->compute_load(part_loaded);
        } else {
            auto *current_cell = dynamic_cast<Cell<Type> *>(*it);

            if ((current_cell->_size.x / (part_loaded->get(POS) - current_cell->_mass_pos).norm()) < BH_THETA)
                current_cell->compute_load(part_loaded);
            else
                update_load(current_cell, part_loaded);
        }
    }
}*/