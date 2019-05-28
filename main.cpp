/**
 * @file main.cpp
 * @author Munier Louis
 * @date 27.05.19
 * @version 1.0
 *
 * Main file of the nBody project.
 */

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 * @include Vector.hpp custom library to have minimal vector implementation
 * @include Tree.hpp library to create a quadtree/octree data structure and interact on different cells/particles
 *
 * @namespace use std to simplify code implementation
 * @namespace use Tree to simplify code implementation
 */
#include "constants.hpp"
#include "Vector.hpp"
#include "Tree.hpp"

using namespace std;
using namespace Tree;

/**
 * Generate NB_PARTICLES particles and call function to store them.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer to initial Cell of the tree to construct
 * @param vec vector to give the Type of the vector to the template function
 * @param list_particles pointer on a list of all the pointers on generated particles (if PRINT defined) to generate csv
 */
template <typename Type>
void generate_data(Cell<Type>* root, Type vec) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5*root->_size.x;
    float y_rnd, p_y = 0.5*root->_size.y;

    vector< Particle<Type>* > other_particle;

    random_device rd;
    uniform_real_distribution<float> dist_m(0, p_m);
    uniform_real_distribution<float> dist_x(-p_x, p_x);
    uniform_real_distribution<float> dist_y(-p_y, p_y);

    for (int i = 0; i < NB_PARTICLES; i++) {
        // Generate random coordinate for a new particle and shift it if it stays in boundaries to stress application
        x_rnd = dist_x(rd);
        if(fabs(x_rnd - SHIFT) < p_x)
            x_rnd -= SHIFT;

        y_rnd = dist_y(rd);
        if(fabs(y_rnd - SHIFT) < p_y)
            y_rnd -= SHIFT;

#if NB_DIM == DIM_3
        float z_rnd, p_z = 0.5*root->_size.z;
        uniform_real_distribution<float> dist_z(-p_z, p_z);

        z_rnd = dist_z(rd);
        if(fabs(z_rnd - SHIFT) < p_z)
            z_rnd -= SHIFT;
#endif

        // Initialize quadtree/octree
        if(root->_prev == root) {
            root->subdivide_tree();
            root->_prev = nullptr;
        }

        // Create new particle
#if NB_DIM == DIM_2
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

        store_particle(root, new_particle, &other_particle);
    }
}

/**
 * Update the load applied to a particle.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param head pointer to the current cell of the tree
 * @param part_loaded pointer to the current particle for which the load is computed
 */
template <typename Type>
void update_load(Cell<Type> *head, Particle<Type> *part_loaded = nullptr) {
    static Cell<Type>* root = head;

    for (auto it = head->_next.begin(); it != head->_next.end(); ++it) {
        if ((*it) == nullptr)
            return;
        else if ((*it)->get_type() == ParticleT) {
            if (part_loaded == nullptr)
                update_load(root, dynamic_cast<Particle<Type> *>(*it));
            else if ((*it) != part_loaded)
                part_loaded->compute_load(dynamic_cast<Particle<Type> *>(*it));//TODO
        } else {
            auto *current_cell = dynamic_cast<Cell<Type> *>(*it);
            if (part_loaded != nullptr && (current_cell->_size.x / (part_loaded->get(POS) - current_cell->_mass_pos).norm()) < BH_THETA)
                current_cell->compute_load(part_loaded); //TODO
            else
                update_load(current_cell, part_loaded);
        }
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
int find_cell_idx(Type origin, Type particle) {
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
 * Store particle in the tree.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param head pointer on current cell to store the particle
 * @param particle pointer on current particle to store
 * @param list_p pointer on a list of pointers on all the particles to store after the current one
 */
template <typename Type>
void store_particle(Cell<Type>* head, Particle<Type>* particle, vector< Particle<Type>* >* list_p) {
    int cell_idx_1 = find_cell_idx(head->_center, particle->get(POS));

    if (head->_next.empty()) {
        head->_next.push_back(particle);
        particle->_parent = head;
        update_cell(particle, true);

        do {
            head = head->_prev;

            if (head == nullptr)
                break;

            update_cell(particle, true);
        } while (head->_prev != nullptr);

        if (!list_p->empty()) {
            particle = list_p->back();
            list_p->pop_back();

            store_particle(head, particle, list_p);
        }
    } else if (head->_next[0]->get_type() == ParticleT) {
        list_p->push_back((Particle<Type>*) head->_next[0]);
        update_cell(dynamic_cast<Particle<Type> *>(head->_next[0]), false);

        /**< clear precedent pointers */
        head->_next.clear();
        particle->_parent = nullptr;

        head->subdivide_tree();
        store_particle(head, particle, list_p);
    } else {
        store_particle((Cell<Type>*) head->_next[cell_idx_1], particle, list_p);
    }
}

/**
 * Update cell total mass and center of mass.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle that need to update its cell
 * @param add boolean value to add or remove properties of the particle
 */
template <typename Type>
void update_cell(Particle<Type> *particle, bool add){
    auto head = dynamic_cast<Cell<Type> *>(particle->_parent);

    if (add) {
        head->_mass_pos = (head->_mass_pos * head->_m + particle->get(POS) * particle->get_mass())
                        / (head->_m + particle->get_mass());
        head->_m += particle->get_mass();
    } else {
        head->_mass_pos = (head->_mass_pos * head->_m - particle->get(POS) * particle->get_mass())
                        / (head->_m - particle->get_mass());
        head->_m -= particle->get_mass();
    }
}

/**
 * Update tree to change place of particles out of original boundaries.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle that change its place
 */
template <typename Type>
void update_tree(Particle<Type> *particle){
    int nb_particles = 0;
    auto parent = dynamic_cast<Cell<Type> *>(particle->_parent);

    vector< Particle<Type>* > other_particle;

    /**< delete pointer from cell to particle */
    parent->_next.clear();

    /**< parent of the particle go back from one level */
    particle->_parent = parent->_prev;
    parent = parent->_prev;

    if (is_out_boundaries(particle)) {
        for (auto it = parent->_next.begin(); it != parent->_next.end(); ++it) {
            if (!dynamic_cast<Cell<Type> *>(*it)->_next.empty())
                nb_particles++;
        }

        /**< delete empty level of the parent node */
        if (nb_particles == 0)
            dynamic_cast<Cell<Type> *>(particle->_parent)->del_level();
    }

    /**< continue to go up in level to find right cell for our particle */
    while(is_out_boundaries(particle)) {
        update_cell(particle, false);
        particle->_parent = dynamic_cast<Cell<Type> *>(particle->_parent)->_prev;

        if (particle->_parent == nullptr) {
            delete particle;
            return;
        }
    }

    store_particle(dynamic_cast<Cell<Type> *>(particle->_parent), particle, &other_particle);
}

/**
 * Test if a particle is out of its cell and return true if it is the case, false otherwise.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle pointer on the tested particle
 * @return true if particle out of its boundaries, false otherwise
 */
template <typename Type>
bool is_out_boundaries(Particle<Type>* particle) {
    Type position = particle->get(POS);
    Type center = dynamic_cast<Cell<Type> *>(particle->_parent)->_center;
    Type cell_size = dynamic_cast<Cell<Type> *>(particle->_parent)->_size;

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
 * If PRINT is defined, generate a csv file to display animation of the result.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param list_particles pointer on a list of pointers on all the particles to write in csv file
 * @param millis_time timestep to change filename and save chronology
 */
#ifdef PRINT
template <typename Type>
void generate_file(Particle<Type>* particle, int millis_time) {
    ofstream csv_file;
    string filename = "../tests/test_" + to_string(millis_time) + ".csv";

    csv_file.open(filename, ios::app);

    /**< Check if file is empty */
    if (csv_file.tellp() == 0) {
#if NB_DIM == DIM_2
        csv_file << "x,y\n";
#elif NB_DIM == DIM_3
        csv_file << "x,y,z\n";
#endif
    }

    csv_file << particle->get(POS).to_file();
    csv_file.close();
}
#endif

template <typename Type>
void update_particles(Cell<Type>* root, int millis_time){
    for (auto it = root->_next.begin(); it != root->_next.end(); ++it) {
        if ((*it) == nullptr) {
            return;
        } else if ((*it)->get_type() == ParticleT) {
            auto current_particle = dynamic_cast<Particle<Type> *>(*it);

            current_particle->update_vel_pos();

            if (is_out_boundaries(current_particle)) {
                update_cell(current_particle, false);
                update_tree(current_particle);
            }

#ifdef PRINT
            generate_file(current_particle, millis_time);
#endif
        } else if ((*it)->get_type() == CellT) {
            update_particles(dynamic_cast<Cell<Type> *>(*it), millis_time);
        }
    }
}

/**
 * Barnes hut part of the algorithm.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param vec_dim vector of the dimensions of the overall area where particles are generated
 */
template <typename Type>
void barnes_hut(Type vec_dim) {
    auto root = new Cell<Type>(Type(), vec_dim, Type());
    root->_prev = root;
    generate_data(root, Type());

    for (int i = 0; i < ITERATIONS; i++) {
        update_load(root);
        update_particles(root, 1000 * i * DELTA_T);
    }
}

/**
 * Main function, compute time to solve problem and store size of the overall area where particle are studied.
 *
 * @param argc default input in c++ main function
 * @param argv default input in c++ main function
 * @return success if no errors are reached
 */
int main(int argc, char *argv[]) {
    auto start = high_resolution_clock::now();
#if NB_DIM == DIM_2
        int width = SIDE, height = SIDE;
        barnes_hut(Vector2f(width, height));
#elif NB_DIM == DIM_3
        int width = SIDE, height = SIDE, depth = SIDE;
        barnes_hut(Vector3f(width, height, depth));
#endif
    auto stop = high_resolution_clock::now();

    cout << duration_cast<microseconds>(stop - start).count() << endl;

    return 0;
}