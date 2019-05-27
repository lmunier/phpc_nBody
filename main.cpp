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

/**
 * If PRINT is defined, generate a csv file to display animation of the result.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param list_particles pointer on a list of pointers on all the particles to write in csv file
 * @param millis_time timestep to change filename and save chronology
 */
#ifdef PRINT
template <typename Type>
void generate_file(vector< Tree::Particle<Type>* > &list_particles, int millis_time) {
    ofstream csv_file;
    string filename = "../tests/test_" + to_string(millis_time) + ".csv";

    csv_file.open(filename);

#if NB_DIM == DIM_2
    csv_file << "x,y\n";
#elif NB_DIM == DIM_3
    csv_file << "x,y,z\n";
#endif

    for (auto it = list_particles.begin(); it != list_particles.end(); ++it) {
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
    vector< Tree::Particle<Type>* > list_particles;

    auto root = new Tree::Cell<Type>(Type(), vec_dim, Type());
    root->_prev = root;
    root->generate_data(list_particles);

    for (int i = 0; i < ITERATIONS; i++) {
        for (auto it = list_particles.begin(); it != list_particles.end(); ++it) {
            //update_load(root, dynamic_cast<Particle<Type> *> (*it));
            dynamic_cast<Tree::Particle<Type> *> (*it)->update_vel_pos();

            /*if (dynamic_cast<Particle<Type> *>(*it)->is_out_boundaries())
                update_tree();
*/
#ifdef PRINT
            generate_file(list_particles, 1000 * i * DELTA_T);
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