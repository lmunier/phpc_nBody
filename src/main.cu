/**
 * @file main.cpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
 *
 * Main file of the nBody project.
 *
 * @mainpage N-Body project
 *
 * Project in the case of the course of Parallel and High Performance Computing at EPFL - Ecole Polytechnique Fédérale
 * de Lausanne.
 *
 * @section Description
 *
 * The n-body problem aims at simulating a dynamical system of particles under the influence of physical forces. We’ll
 * restrain on the gravity field applied on celestial bodies :
 * @f[ F_{ij} = \frac{G m_i m_j (q_j - q_i)}{|| q_j - q_i ||} @f]
 *
 * where @f$ G @f$ is the gravitational constant, @f$ m_i @f$ and @f$ m_j @f$ the masses of the @f$ i @f$-th and
 * @f$j@f$-th bodies and @f$ q_i @f$ and @f$ q_j @f$ their positions.
 *
 * @section Implementation
 *
 * The solution is implemented using Barnes-Hut algorithm and quadtree/octree data-structure. The particularity is that
 * the user can choose differents parameters as define :
 * - NB_DIM : set at DIM_2 or DIM_3 to choose between 2D and 3D problem solution
 * - NB_PARTICLES : number of particles
 * - PRINT : save solution in csv file, at each time-step, to display solution animation in programs (e.g. paraview)
 * - DELTA_T : time-step of each iteration
 */

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 * @include Vector.hpp custom library to have minimal vector implementation
 * @include Tree.hpp library to create a quadtree/octree data structure and interact on different cells/particles
 */
#include "constants.hpp"
#include "Cell.hpp"
#include "Particle.hpp"

/**
 * @namespace std to simplify code implementation
 * @namespace Tree to simplify code implementation
 */
using namespace std;
using namespace Tree;

/**
 * Generate NB_PARTICLES particles and call function to store them.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer to initial Cell of the tree to construct
 * @param vec vector to give the Type of the vector to the template function
 */
template <typename Type>
void generate_data(Cell<Type>* root, Type vec) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5f * OCCUPATION_PERC * root->get(DIM).x;
    float y_rnd, p_y = 0.5f * OCCUPATION_PERC * root->get(DIM).y;

    vector< AbstractType<Type>* > other_particle{};

    random_device rd;
    uniform_real_distribution<float> dist_m(0, p_m);
    uniform_real_distribution<float> dist_x(-p_x, p_x);
    uniform_real_distribution<float> dist_y(-p_y, p_y);

    for (int i = 0; i < NB_PARTICLES; i++) {
        /** Generate random coordinate for a new particle and shift it if it stays in boundaries
         * to stress application */
        x_rnd = dist_x(rd);
        if(fabs(x_rnd - SHIFT) < p_x)
            x_rnd -= SHIFT;

        y_rnd = dist_y(rd);
        if(fabs(y_rnd - SHIFT) < p_y)
            y_rnd -= SHIFT;

#if NB_DIM == DIM_3
        float z_rnd, p_z = 0.5f * OCCUPATION_PERC * root->get(DIM).z;
        uniform_real_distribution<float> dist_z(-p_z, p_z);

        z_rnd = dist_z(rd);
        if(fabs(z_rnd - SHIFT) < p_z)
            z_rnd -= SHIFT;
#endif

        /** Initialize quadtree/octree */
        if(root->get_parent() == root) {
            root->subdivide_tree();
            root->set_parent(nullptr);
        }

        /** Create new particle */
#if NB_DIM == DIM_2
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

        root->store_particle(new_particle, nullptr);
    }
}

/**
 * Update the load applied to a particle by implementation a Depth First Search on the quadtree/octree data-structure.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param head pointer to the current cell of the tree
 * @param part_loaded pointer to the current particle for which the load is computed
 */
template <typename Type>
void update_load(AbstractType<Type> *head, AbstractType<Type> *part_loaded = nullptr) {
    static AbstractType<Type>* root = head;

    for (auto n : head->get_next()) {
        /** If next element is empty */
        if (n == nullptr) {
            return;
        } /** If next element is a particle */
        else if (n->get_type() == ParticleT) {
            if (part_loaded == nullptr)
                update_load(root, n);
            else if (n != part_loaded)
                n->compute_load(part_loaded);
        } /** If next element is a cell */
        else {
            if (part_loaded != nullptr && (n->get(DIM).x / (part_loaded->get(POS) - n->get(MASS_POS)).norm()) < BH_THETA)
                n->compute_load(part_loaded);
            else
                update_load(n, part_loaded);
        }
    }
}

/**
 * Update position and velocity for each particle. Generate a csv file with all the position if needed.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer on the node of the previous cell
 * @param iter current iteration of the solution
 */
template <typename Type>
void update_particles_pos(AbstractType<Type>* root, int iter, const string& dir){
    for (auto n : root->get_next()) {
        /** If next element is empty */
        if (n == nullptr) {
            return;
        } else if (n->get_type() == ParticleT) {
            n->update_vel_pos();

#ifdef PRINT
            if (n != nullptr)
                generate_file(n, 1000 * iter * DELTA_T, dir);
#endif
        } else if (n->get_type() == CellT) {
            update_particles_pos(n, iter, dir);
        }
    }
}

/**
 * Update the current quadtree/octree datastructure with the new position of each particle if they pass to an other
 * cell.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer on the node of the previous cell
 */
template <typename Type>
void update_particles_tree(AbstractType<Type>* root){
    auto next = root->get_next();

    for (auto it = next.begin(); it != next.end(); ++it) {
        /** If next element is empty or is already deleted */
        if ((*it) == nullptr || !(*it)->get_state()) {
            return;
        } /** If next element is a particle */
        else if ((*it)->get_type() == ParticleT) {
            if ((*it)->is_out_boundaries()) {
                if ((*it)->update_tree() == -1)
                    delete *it;
            }
        } /** If next element is a cell */
        else if ((*it)->get_type() == CellT) {
            update_particles_tree(*it);
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
void barnes_hut(Type vec_dim, const string& dir) {
    auto root = new Cell<Type>(Type(), vec_dim, Type());
    root->set_parent(root);
    generate_data(root, Type());

    for (int i = 1; i <= ITERATIONS; i++) {
        update_load(root);
        update_particles_pos(root, i, dir);
        update_particles_tree(root);
    }
}

/**
 * If PRINT is defined, generate a csv file to display animation of the result in external software (e.g. paraview).
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle pointer on the particle to write in csv file
 * @param millis_time timestep to change filename and save chronology
 */
#ifdef PRINT
template <typename Type>
void generate_file(AbstractType<Type>* particle, int millis_time, const string& dir) {
    ofstream csv_file;
    string filename = dir + "/out_" + to_string(millis_time) + ".csv";

    csv_file.open(filename, ios::app);

    /** Check if file is empty to write the title of each column */
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

/**
 * Main function, compute time to solve problem and store size of the overall area where particle are studied.
 *
 * @param argc default input in c++ main function
 * @param argv default input in c++ main function
 * @return success if no errors are reached
 */
int main(int argc, char *argv[]) {
    /** Find/set the device. */
    /** The test requires an architecture SM35 or greater (CDP capable). */
    int cuda_device = findCudaDevice(argc, (const char **)argv);
    cudaDeviceProp deviceProps;
    checkCudaErrors(cudaGetDeviceProperties(&deviceProps, cuda_device));
    int cdpCapable = (deviceProps.major == 3 && deviceProps.minor >= 5) || deviceProps.major >=4;

    printf("GPU device %s has compute capabilities (SM %d.%d)\n", deviceProps.name, deviceProps.major, deviceProps.minor);

    if (!cdpCapable)
    {
        std::cerr << "cdpQuadTree requires SM 3.5 or higher to use CUDA Dynamic Parallelism.  Exiting...\n" << std::endl;
        exit(EXIT_SUCCESS);
    }

    int width = SIDE, height = SIDE;
    auto start = high_resolution_clock::now();

#if NB_DIM == DIM_2
    if (argv[1])
        barnes_hut(Vector2f(width, height), argv[1]);
    else
        barnes_hut(Vector2f(width, height), "../output");
#elif NB_DIM == DIM_3
    int depth = SIDE;

    if (argv[1])
        barnes_hut(Vector3f(width, height, depth), argv[1]);
    else
        barnes_hut(Vector3f(width, height, depth), "../output");
#endif
    auto stop = high_resolution_clock::now();

    cout << duration_cast<microseconds>(stop - start).count() << endl;

    return 0;
}