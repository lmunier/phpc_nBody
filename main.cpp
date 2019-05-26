#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "constants.hpp"
#include "Vector.hpp"
#include "Tree.hpp"

using namespace std;
using namespace Tree;

template <typename Type>
void generate_data(Cell<Type>* root, Type vec, int nb_particles, vector< Particle<Type>* >* list_particles = nullptr) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5*root->_size.x;
    float y_rnd, p_y = 0.5*root->_size.y;

    vector< Particle<Type>* > other_particle;

    random_device rd;
    uniform_real_distribution<float> dist_m(0, p_m);
    uniform_real_distribution<float> dist_x(-p_x, p_x);
    uniform_real_distribution<float> dist_y(-p_y, p_y);

    for (int i = 0; i < nb_particles; i++) {
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
            subdivide_tree(root, vec);
            root->_prev = nullptr;
        }

        // Create new particle
#if NB_DIM == DIM_2
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

#ifdef PRINT
        list_particles->push_back(new_particle);
#endif

        store_particle(root, new_particle, &other_particle);
    }
}

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

template <typename Type>
void update_vel_pos(Particle<Type> *part_loaded) {
    Type new_velocity = part_loaded->get(LOAD) * (DELTA_T/part_loaded->get_mass()) + part_loaded->get(VEL);
    Type new_position = new_velocity * DELTA_T + part_loaded->get(POS);

    part_loaded->set(VEL, new_velocity);
    part_loaded->set(POS, new_position);
    part_loaded->set(LOAD, Type());
}

template <typename Type>
void subdivide_tree(Cell<Type>* previous, Type vec) {
    bool y = true, z = false;
    Type size = previous->_size*0.5;

    for(int n = 0; n < pow(2, NB_DIM); n++) {
        if (n%2 == 0)
            y = !y;

        if (n == 4)
            z = !z;

#if NB_DIM == DIM_2
        Type center = Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)));
#elif NB_DIM == DIM_3
        Type center = Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)), size.z*(0.5*pow(-1, z)));
#endif

        center = previous->_center + center;

        auto next = new Cell<Type>(center, size, center);
        previous->_next.push_back(next);
        next->_prev = previous;
    }
}

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

template <typename Type>
void store_particle(Cell<Type>* head, Particle<Type>* particle, vector< Particle<Type>* >* list_p = nullptr) {
    int cell_idx_1 = find_cell_idx(head->_center, particle->get(POS));

    if (head->_next.empty()) {
        head->_next.push_back(particle);
        head->_m += particle->get_mass();
        head->_mass_pos = head->_mass_pos + particle->get(POS);

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

        subdivide_tree(head, Type());
        store_particle(head, particle, list_p);
    } else {
        store_particle((Cell<Type>*) head->_next[cell_idx_1], particle, list_p);
    }
}

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

template <typename Type>
void barnes_hut(Type vec_dim) {
    vector< Particle<Type>* > list_particles;

    auto root = new Cell<Type>(Type(), vec_dim, Type());
    root->_prev = root;
    generate_data(root, Type(), NB_PARTICLES, &list_particles);

    for (int i = 0; i < ITERATIONS; i++) {
        for (auto it = list_particles.begin(); it != list_particles.end(); ++it) {
            update_load(root, *it);
            update_vel_pos(*it);

#ifdef PRINT
            generate_file(&list_particles, 1000 * i * DELTA_T);
#endif
        }
    }
}

int main(int argc, char *argv[]) {
#if NB_DIM == DIM_2
        int width = SIDE, height = SIDE;
        barnes_hut(Vector2f(width, height));
#elif NB_DIM == DIM_3
        int width = SIDE, height = SIDE, depth = SIDE;
        barnes_hut(Vector3f(width, height, depth));
#endif

    return 0;
}