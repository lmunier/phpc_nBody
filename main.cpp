#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <type_traits>
#include <typeinfo>

#include "constants.hpp"
#include "Vector.hpp"
#include "Tree.hpp"

using namespace std;
using namespace Tree;

template <typename Type>
Type compute_load(Particle<Type>* part1, Particle<Type>* part2) {
    Type tmp = part1->get(POS) - part2->get(POS);
    return tmp*(G*part1->get_mass()*part2->get_mass()/tmp.norm());
}

template <typename Type>
void generate_data(Cell<Type>* root, Type vec, int nb_particles) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5*root->_size.x;
    float y_rnd, p_y = 0.5*root->_size.y;

    random_device rd;
    uniform_real_distribution<float> dist_m(0, p_m);
    uniform_real_distribution<float> dist_x(-p_x, p_x);
    uniform_real_distribution<float> dist_y(-p_y, p_y);

    for (int i = 0; i < nb_particles; i++) {
        // Generate random coordinate for a new particle and shift it if it stays in boundaries to stress application
        x_rnd = dist_x(rd);
        if(fabs(dist_x(rd) - SHIFT) < p_x)
            x_rnd -= SHIFT;

        y_rnd = dist_y(rd);
        if(fabs(dist_y(rd) - SHIFT) < p_y)
            y_rnd -= SHIFT;

#if NB_DIM == DIM_3
        float z_rnd, p_z = 0.5*root->_size.z;
        uniform_real_distribution<float> dist_z(-p_z, p_z);
        z_rnd = dist_z(rd);
#endif

        // Initialize quadtree/octree
        if(root->_prev == root) {
            subdivide_tree(root, vec);
            root->_prev = nullptr;
        }

        // Create new particle
#if NB_DIM == DIM_2
        auto new_particle = new Particle<vec_Type>(dist_m(rd), vec_Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

        store_particle(root, new_particle);

        // Test
        new_particle->get(POS).print();
    }
}

template <typename Type>
void subdivide_tree(Cell<Type>* previous, Type vec) {
    bool y = true, z = false;

    for(int n = 0; n < pow(2, NB_DIM); n++) {
        Type size = previous->_size*0.5;

        if (n%2 == 0)
            y = !y;

        if (n == 4)
            z = !z;

#if NB_DIM == DIM_2
        vec_Type center = vec_Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)));
#elif NB_DIM == DIM_3
        Type center = Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)), size.z*(0.5*pow(-1, z)));
#endif

        auto next = new  Cell<Type>(center, size, center);
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
void store_particle(Cell<Type>* head, Particle<Type>* particle_1, Particle<Type>* particle_2 = nullptr) {
    int cell_idx_1 = find_cell_idx(head->_center, particle_1->get(POS));

    if (particle_2 != nullptr)
        int cell_idx_2 = find_cell_idx(head->_center, particle_2->get(POS));

    if (head->_next.empty())
        head->_next.push_back(particle_1);
    else if (head->_next[cell_idx_1]->get_type() == CellT)
        store_particle((Cell<Type>*) head->_next[cell_idx_1], particle_1, particle_2);
    else {
        particle_2 = (Particle<Type>*) head->_next[cell_idx_1];
        subdivide_tree((Cell<Type>*) head->_next[cell_idx_1], Type());
        store_particle(head, particle_1, particle_2);
    }
}


template <typename Type>
void barnes_hut(Type vec_dim) {
    auto root = new Cell<Type>(Type(), vec_dim, Type());
    root->_prev = root;
    generate_data(root, Type(), NB_PARTICLES);

    // Test vector implementation
    Particle<Type> ptr = Particle<Type>(2.0f, Type(2.0f, 5.0f));
    Particle<Type> ptr2 = Particle<Type>(2.0f, Type(3.0f, -6.0f));
    ptr2.get(POS).print();

    Particle<Type> ptr3 = Particle<Type>();
    ptr3.set(POS, ptr.get(POS) + ptr2.get(POS));
    ptr3.get(POS).print();

    // Compute load
    compute_load(&ptr, &ptr2).print();
}

int main(int argc, char *argv[]) {
#if NB_DIM == DIM_2
        int width = 1920, height = 1200;
        barnes_hut(Vector2f(width, height));
#elif NB_DIM == DIM_3
        int width = 1920, height = 1200, depth = 1200;
        barnes_hut(Vector3f(width, height, depth));
#endif

    return 0;
}