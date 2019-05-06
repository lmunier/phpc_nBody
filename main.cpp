#include <iostream>
#include "Particle.hpp"

int main(int argc, char *argv[]) {
    int NB_DIM = DIM_2;

    if(argc != 0)
        NB_DIM = (unsigned int) *argv[0];

    if (NB_DIM == DIM_2) {
        auto ptr = new Particle<Vector2f>(2.0f, Vector2f(2.0f, 5.0f));
        ptr->get_pos().print();
    } else {
        auto ptr = new Particle<Vector3f>(2.0f, Vector3f(2.0f, 5.0f, 4.0f));
        ptr->get_pos().print();
    }

    return 0;
}