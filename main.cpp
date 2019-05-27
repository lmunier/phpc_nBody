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
void generate_data(Cell<Type>* root, Type vec, vector< Particle<Type>* >* list_particles = nullptr) {
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
        if(i == 0)
            new_particle->set(MASS, 1000);

#ifdef PRINT
        list_particles->push_back(new_particle);
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
void update_load(Cell<Type> *head, Particle<Type> *part_loaded) {
    for (auto it = head->_next.begin(); it != head->_next.end(); ++it) {
        if ((*it) == nullptr)
            return;
        else if ((*it)->get_type() == ParticleT) {
            if ((*it) != part_loaded)
                part_loaded->compute_load(dynamic_cast<Particle<Type> *>(*it));
        } else
            update_load(dynamic_cast<Cell<Type> *> (*it), part_loaded);
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
void store_particle(Cell<Type>* head, Particle<Type>* particle, vector< Particle<Type>* >* list_p = nullptr) {
    int cell_idx_1 = find_cell_idx(head->_center, particle->get(POS));

    if (head->_next.empty()) {
        head->_next.push_back(particle);
        head->_m += particle->get_mass();
        head->_mass_pos = (head->_mass_pos * head->_m + particle->get(POS) * particle->get_mass()) *
                          ( 1 / (head->_m + particle->get_mass()));

        do {
            head = head->_prev;

            head->_m += particle->get_mass();
            head->_mass_pos = head->_mass_pos + particle->get(POS);
        } while (head->_prev != nullptr);

        if (!list_p->empty()) {
            particle = list_p->back();
            list_p->pop_back();

            store_particle(head, particle, list_p);
        }
    } else if (head->_next[0]->get_type() == ParticleT) {
        list_p->push_back((Particle<Type>*) head->_next[0]);
        head->_next.clear();

        head->subdivide_tree();
        store_particle(head, particle, list_p);
    } else {
        store_particle((Cell<Type>*) head->_next[cell_idx_1], particle, list_p);
    }
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
void generate_file(vector< Particle<Type>* >* list_particles, int millis_time) {
    ofstream csv_file;
    string filename = "../tests/test_" + to_string(millis_time) + ".csv";

    csv_file.open(filename);

#if NB_DIM == DIM_2
    csv_file << "x,y\n";
#elif NB_DIM == DIM_3
    csv_file << "x,y,z\n";
#endif

    for (auto it = list_particles->begin(); it != list_particles->end(); ++it) {
        csv_file << (*it)->get(POS).to_file();
    }

    csv_file.close();
}
#endif

/**
 * Barnes hut part of the algorithm.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param vec_dim vector of the dimensions of the overall area where particles are generated
 */
template <typename Type>
void barnes_hut(Type vec_dim) {
    vector< Particle<Type>* > list_particles;

    auto root = new Cell<Type>(Type(), vec_dim, Type());
    root->_prev = root;
    generate_data(root, Type(), &list_particles);

    for (int i = 0; i < ITERATIONS; i++) {
        for (auto it = list_particles.begin(); it != list_particles.end(); ++it) {
            update_load(root, *it);
            dynamic_cast<Particle<Type> *> (*it)->update_vel_pos();

#ifdef PRINT
            generate_file(&list_particles, 1000 * i * DELTA_T);
#endif
        }
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