#include <iostream>
#include "Particle.hpp"
#include "constants.hpp"


int main(int argc, char *argv[]) {
    int NB_DIM = DIM_2;

    if(argc == 2)
        NB_DIM = std::stoi(argv[1]);

    std::cout << NB_DIM << std::endl;

    if (NB_DIM == DIM_2) {
        Particle<Vector2f> ptr = Particle<Vector2f>(2.0f, Vector2f(2.0f, 5.0f));
        Particle<Vector2f> ptr2 = Particle<Vector2f>(2.0f, Vector2f(2.0f, -6.0f));
        ptr.get(POS).print_pos();

        Particle<Vector2f> ptr3 = Particle<Vector2f>();
        ptr3.set(POS, ptr.get(POS) + ptr2.get(POS));
        ptr3.get(POS).print_pos();
    } else {
        Particle<Vector3f> ptr = Particle<Vector3f>(2.0f, Vector3f(2.0f, 5.0f, 9.0f));
        Particle<Vector3f> ptr2 = Particle<Vector3f>(2.0f, Vector3f(2.0f, 6.0f, -7.0f));
        ptr.get(POS).print_pos();

        Particle<Vector2f> ptr3 = Particle<Vector2f>();
        ptr3.set(POS, ptr.get(POS) + ptr2.get(POS));
        ptr3.get(POS).print_pos();
    }

    return 0;
}