#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <type_traits>

#include "constants.hpp"
#include "Vector.hpp"
#include "Particle.hpp"
#include "Cell.hpp"

using namespace std;

template <typename Type>
Type compute_load(Particle<Type>* part1, Particle<Type>* part2) {
    Type tmp = part1->get(POS) - part2->get(POS);
    return tmp*(G*part1->get_mass()*part2->get_mass()/tmp.norm());
}

template <typename cell_Type, typename vec_Type>
void generate_data(cell_Type root, vec_Type vec, int nb_particles) {
    int x_rnd, p_x = root->_size.x;
    int y_rnd, p_y = root->_size.y;

    random_device rd;
    uniform_int_distribution<int> dist_x(-p_x, p_x);
    uniform_int_distribution<int> dist_y(-p_y, p_y);

    for (int i = 0; i < nb_particles; i++) {
        x_rnd = dist_x(rd) - SHIFT;
        y_rnd = dist_y(rd) - SHIFT;

#if NB_DIM == DIM_3
        if(is_same<Vector3f, vec_Type>::value) {
            int z_rnd, p_z = root->_size.z;
            uniform_int_distribution<int> dist_z(-p_z, p_z);
            z_rnd = dist_z(rd);
        }
#endif
        if(root->_prev == nullptr) {
            bool y = true, z = false;

            for(int n = 0; n < pow(2, NB_DIM); n++) {
                vec_Type size = root->_size*0.5;

                if (n%2 == 0)
                    y = !y;

                if (n == 4)
                    z = !z;
#if NB_DIM == DIM_2
                vec_Type m_pos = vec_Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)));
                root->_next.push_back(new Cell<vec_Type>(size, m_pos));
#elif NB_DIM == DIM_3
                vec_Type m_pos = vec_Type(size.x*(0.5*pow(-1, (n+1)%2)), size.y*(0.5*pow(-1, y)), size.z*(0.5*pow(-1, z)));
                root->_next.push_back(new Cell<vec_Type>(size, m_pos));
#endif
                m_pos.print();
            }
        }
    }
}

template <typename Type>
void barnes_hut(Type vec_dim) {
    auto root = new Cell<Type>(vec_dim, Type());
    root->_prev = nullptr;
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