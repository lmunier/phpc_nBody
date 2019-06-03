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
void generate_data(thrust::device_vector< Particle<Type>* >* particles, Type size) {
    float p_m = MASS_MAX;
    float x_rnd, p_x = 0.5f * OCCUPATION_PERC * size.x;
    float y_rnd, p_y = 0.5f * OCCUPATION_PERC * size.y;

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
        float z_rnd, p_z = 0.5f * OCCUPATION_PERC * size.z;
        uniform_real_distribution<float> dist_z(-p_z, p_z);

        z_rnd = dist_z(rd);
        if(fabs(z_rnd - SHIFT) < p_z)
            z_rnd -= SHIFT;
#endif

        /** Create new particle */
#if NB_DIM == DIM_2
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd));
#elif NB_DIM == DIM_3
        auto new_particle = new Particle<Type>(dist_m(rd), Type(x_rnd, y_rnd, z_rnd));
#endif

        particles->push_back(new_particle);
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
void update_load(thrust::device_vector< Particle<Type>* >* particles) {
    for (auto p_i : *particles) {
        Type load = Type();

        for (auto p_j : *particles) {
            Type tmp = (*p_j).get(POS) - (*p_i).get(POS);
            float d = max(tmp.norm(), EPSILON);
            load = load + tmp * (G * (*p_i).get_mass() * (*p_j).get_mass()) / d;
        }

        (*p_i).set(LOAD, load);
    }
}

struct square {
    __host__ __device__ float operator()(float x) {
        return x * x;
    }
};
    
__host__ __device__ float snrm2_fast (thrust::device_vector<float>& x) {
    // with fusion
    return sqrt( thrust::transform_reduce(x.begin(), x.end(), square(), 0.0f, thrust::plus<float>()));
}

/**
 * Update position and velocity for each particle. Generate a csv file with all the position if needed.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer on the node of the previous cell
 * @param iter current iteration of the solution
 */
template <typename Type>
void update_particles_pos(vector< Particle<Type>* >* particles, Type dim, int iter, const string& dir){
    for (auto p : *particles) {
        p->update_vel_pos();

        if (p->is_out_boundaries(dim))
            delete p;

#ifdef PRINT
        generate_file(p, 1000 * iter * DELTA_T, dir);
#endif
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
    /*int width = SIDE, height = SIDE;
    auto start = high_resolution_clock::now();

    string dir = "";

    if (argv[1])
        dir = argv[1];
    else
        dir = "../output";*/

    thrust::device_vector<float> d_vec(10);
    thrust::generate(d_vec.begin(), d_vec.end(), rand);

    float result = snrm2_fast(d_vec);
    cudaDeviceSynchronize();

    //printf("%f\n", result);

/*#if NB_DIM == DIM_2
    Vector2f size = Vector2f(width, height);
    thrust::device_vector< Particle<Vector2f>* > particles{};
    generate_data(&particles, size);

    for (int i = 1; i <= ITERATIONS; i++) {
        update_load(&particles);
        //update_particles_pos(&particles, size, i, dir);
    }
#elif NB_DIM == DIM_3
    int depth = SIDE;

    Vector3f size = Vector3f(width, height, depth);
    thrust::device_vector< Particle<Vector3f>* > particles{};
    generate_data(&particles, Vector3f(width, height, depth));

    for (int i = 1; i <= ITERATIONS; i++) {
        update_load(&particles);
        //update_particles_pos(&particles, size, i, dir);
    }
#endif*/
    //auto stop = high_resolution_clock::now();

    /** Print all the parameters */
    /*cout << "Brut force" << endl;
    cout << "Epsilon " << EPSILON << endl;
    cout << "Nb particles " << NB_PARTICLES << endl;
    cout << "Nb dimensions " << NB_DIM << endl;
    cout << "Side " << SIDE << endl;
    cout << "Shift " << SHIFT << endl;
    cout << "Occupation percentage " << OCCUPATION_PERC << endl;
    cout << "Maximum mass " << MASS_MAX << endl;
    cout << "Delta t " << DELTA_T << endl;
    cout << "Nb iterations " << ITERATIONS << endl;

    cout << duration_cast<microseconds>(stop - start).count() << endl;*/

    return 0;
}