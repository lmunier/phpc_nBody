//
// Created by lmunier on 12.05.19.
//

#ifndef PROJECT_CELL_HPP
#define PROJECT_CELL_HPP

#include <vector>

#include "constants.hpp"
#include "Particle.hpp"
#include "Vector.hpp"

using namespace std;

template<typename Type>
class Cell : public Particle<Type> {
public:
    explicit Cell(Type dim = Type(), Type mass_pos = Type()) {
        this->_size = dim;
        this->_mass_pos = mass_pos;
    }

    void del_element(Particle<Type>* p) {
        p->del();

        for(auto it = this->_next.begin(); it != this->_next.end(); ++it) {
            if (it == p) {
                (*it) = nullptr;
                break;
            }
        }
    }

    Type _size;
    Type _mass_pos;

    Cell* _prev;
    vector<Particle<Type>*> _next{};
};


#endif //PROJECT_CELL_HPP
